from adam_community.tool import Tool
from datetime import datetime
from dataclasses import dataclass, asdict
from pathlib import Path
import gzip
import json
import os
import re
import shutil
import subprocess
import tempfile
from typing import Any, Dict, List, Optional, Tuple


@dataclass
class QCConfig:
    """QC 配置对象"""
    assay: str = "RNA"
    qc_method: str = "manual"
    mt_pattern: str = "^MT-"
    mt_col_name: str = "percent.mt"
    min_features: int = 200
    max_features: int = 2500
    min_counts: int = 0
    max_counts: int = 0
    max_percent_mt: int = 5
    save_plots: bool = True
    plot_format: str = "both"
    output_object: str = "seurat"


class QCError(Exception):
    """QC 相关错误"""
    pass


class InputResolver:
    """输入解析器：读取 standardized_object.rds"""

    def __init__(self, work_dir: Path):
        self.work_dir = Path(work_dir)
        self.upload_dir = self.work_dir / "uploaded_raw"
        self.upload_dir.mkdir(parents=True, exist_ok=True)

    def resolve(self, kwargs: Dict[str, Any], config: QCConfig) -> Dict[str, Any]:
        """解析输入文件"""

        # 1. 优先从网页表单/上传对象解析
        input_path = self._resolve_uploaded_rds(kwargs)

        # 2. 兜底：本地调试时检查当前目录 standardized_object.rds
        if input_path is None:
            direct_path = Path("standardized_object.rds")
            if direct_path.exists() and direct_path.is_file():
                input_path = direct_path.resolve()

        # 3. 再兜底：当前目录任意 .rds
        if input_path is None:
            cwd = Path.cwd()
            rds_files = list(cwd.glob("*.rds"))
            if rds_files:
                for f in rds_files:
                    if f.name == "standardized_object.rds":
                        input_path = f.resolve()
                        break
                if input_path is None:
                    input_path = rds_files[0].resolve()

        if input_path is None:
            raise QCError("未找到输入对象文件 standardized_object.rds（或其他 .rds 文件）")

        if not input_path.exists():
            raise QCError(f"输入文件不存在: {input_path}")

        if not input_path.is_file():
            raise QCError(f"输入路径不是文件: {input_path}")

        if input_path.suffix.lower() != ".rds":
            raise QCError(f"输入文件不是 .rds 格式: {input_path.name}")

        return {
            "input_file": str(input_path),
            "input_path": str(input_path),
            "file_exists": True,
            "file_size": input_path.stat().st_size,
        }

    def _resolve_uploaded_rds(self, kwargs: Dict[str, Any]) -> Optional[Path]:
        """优先解析网页端上传的 .rds 文件"""
        candidates = self._collect_uploaded_paths(kwargs)

        if not candidates:
            return None

        # 优先 standardized_object.rds
        for p in candidates:
            if p.name == "standardized_object.rds" and p.suffix.lower() == ".rds":
                return p

        # 其次任意 .rds
        for p in candidates:
            if p.suffix.lower() == ".rds":
                return p

        return None

    def _collect_uploaded_paths(self, kwargs: Dict[str, Any]) -> List[Path]:
        """收集上传文件路径，优先处理 args['input_file']，兼容 kwargs['files']"""
        paths: List[Path] = []
        args = kwargs.get("args", {}) or {}

        # 1) 处理网页表单 input_file
        value = args.get("input_file")
        if value is not None:
            parsed = self._parse_uploaded_item(value, default_name="standardized_object.rds")
            if parsed is not None:
                paths.append(parsed)

        # 2) 向后兼容 kwargs["files"]
        for item in kwargs.get("files", []) or []:
            parsed = self._parse_uploaded_item(item, default_name="uploaded.rds")
            if parsed is not None:
                paths.append(parsed)

        # 去重
        deduped: List[Path] = []
        seen = set()
        for p in paths:
            sp = str(p.resolve()) if p.exists() else str(p)
            if sp not in seen:
                seen.add(sp)
                deduped.append(p)

        return deduped

    def _parse_uploaded_item(self, item: Any, default_name: str) -> Optional[Path]:
        """解析单个上传对象"""

        # 情况 1：直接是字符串路径
        if isinstance(item, str) and item.strip():
            path = Path(item)
            if path.exists() and path.is_file():
                return path.resolve()
            return None

        # 情况 2：字典对象（网页端常见）
        if isinstance(item, dict):
            # 先尝试真实路径字段
            for key in ("path", "file_path", "filepath", "local_path", "saved_path", "file"):
                candidate = item.get(key)
                if isinstance(candidate, str) and candidate.strip():
                    path = Path(candidate)
                    if path.exists() and path.is_file():
                        return path.resolve()

            # 如果只有 content，则落盘
            if "content" in item:
                original_name = (
                    item.get("filename")
                    or item.get("name")
                    or default_name
                )
                return self._save_uploaded_content(item["content"], str(original_name))

            # 注意：filename 只当文件名，不直接当本地路径
            filename = item.get("filename") or item.get("name")
            if isinstance(filename, str) and filename.strip():
                candidate = Path(filename)
                if candidate.exists() and candidate.is_file():
                    return candidate.resolve()

        return None

    def _save_uploaded_content(self, content: Any, original_name: str) -> Path:
        """将网页端传入的内容保存为本地文件"""
        import base64
        import binascii
        import os
        import re

        safe_name = os.path.basename(original_name) if original_name else "uploaded.rds"
        safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", safe_name)

        if not safe_name.lower().endswith(".rds"):
            safe_name += ".rds"

        out_path = self.upload_dir / safe_name

        if isinstance(content, str):
            try:
                file_bytes = base64.b64decode(content, validate=True)
            except (binascii.Error, ValueError):
                file_bytes = content.encode("utf-8")
        elif isinstance(content, bytes):
            file_bytes = content
        else:
            file_bytes = str(content).encode("utf-8")

        with open(out_path, "wb") as f:
            f.write(file_bytes)

        return out_path.resolve()

class RPipelineBuilder:
    """构建并运行 R 脚本进行 QC 分析"""
    
    def __init__(self, work_dir: Path):
        self.work_dir = Path(work_dir)
    
    def write_script(self) -> Path:
        """写入 R 脚本"""
        script_path = self.work_dir / "run_qc.R"
        script_path.write_text(self._script_text(), encoding="utf-8")
        return script_path
    
    def _script_text(self) -> str:
        """R 脚本内容"""
        return r'''args <- commandArgs(trailingOnly = TRUE)
config_path <- args[[1]]
input_path <- args[[2]]
output_dir <- args[[3]]

read_json <- function(path) {
  txt <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  jsonlite::fromJSON(txt, simplifyVector = FALSE)
}

cfg <- read_json(config_path)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# 1. 读取输入对象
message("[R] 读取输入对象: ", input_path)
obj <- readRDS(input_path)

# 检查对象类型
if (!inherits(obj, "Seurat")) {
  stop("当前版本仅支持 Seurat 对象。输入对象类型: ", class(obj)[1])
}

# 检查 assay
if (!cfg$assay %in% names(obj@assays)) {
  stop("指定的 assay '", cfg$assay, "' 在对象中不存在。可用的 assay: ", paste(names(obj@assays), collapse = ", "))
}

DefaultAssay(obj) <- cfg$assay

# 2. 准备 metadata 和计算 QC 指标
message("[R] 计算 QC 指标...")

# 检查并计算 nFeature_RNA
if (!"nFeature_RNA" %in% colnames(obj[[]])) {
  message("[R] nFeature_RNA 不存在，从 counts 矩阵重新计算")
  obj$nFeature_RNA <- colSums(GetAssayData(obj, assay = cfg$assay, slot = "counts") > 0)
} else {
  message("[R] 使用现有的 nFeature_RNA")
}

# 检查并计算 nCount_RNA
if (!"nCount_RNA" %in% colnames(obj[[]])) {
  message("[R] nCount_RNA 不存在，从 counts 矩阵重新计算")
  obj$nCount_RNA <- colSums(GetAssayData(obj, assay = cfg$assay, slot = "counts"))
} else {
  message("[R] 使用现有的 nCount_RNA")
}

# 计算线粒体比例
mt_pattern <- cfg$mt_pattern
mt_col_name <- cfg$mt_col_name

if (!mt_col_name %in% colnames(obj[[]])) {
  message("[R] 计算线粒体比例: pattern = ", mt_pattern)
  mt_genes <- grep(mt_pattern, rownames(obj), value = TRUE, ignore.case = TRUE)
  
  if (length(mt_genes) == 0) {
    warning("未找到匹配线粒体模式 '", mt_pattern, "' 的基因。percent.mt 将设为 0。")
    obj[[mt_col_name]] <- 0
  } else {
    message("[R] 找到 ", length(mt_genes), " 个线粒体基因")
    obj[[mt_col_name]] <- PercentageFeatureSet(obj, pattern = mt_pattern, assay = cfg$assay)
  }
} else {
  message("[R] 使用现有的 ", mt_col_name)
}

# 3. 生成 QC 前指标表
message("[R] 生成 QC 前指标表...")
qc_metrics <- data.frame(
  cell_id = colnames(obj),
  nFeature_RNA = obj$nFeature_RNA,
  nCount_RNA = obj$nCount_RNA,
  percent_mt = obj[[mt_col_name, drop = TRUE]],
  pass_qc = NA,
  remove_reason = "",
  stringsAsFactors = FALSE
)

write.csv(qc_metrics, gzfile(file.path(output_dir, "qc_metrics_all_cells.csv.gz")), row.names = FALSE)

# 4. 生成 QC 前图
if (isTRUE(cfg$save_plots)) {
  message("[R] 生成 QC 前图...")
  plots_dir <- file.path(output_dir, "plots")
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Violin plot
  vln_before <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", mt_col_name), 
                        pt.size = 0.1, ncol = 3)
  
  # Scatter plots
  scatter1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = mt_col_name)
  scatter2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  # 保存图像
  plot_format <- cfg$plot_format
  if (plot_format %in% c("pdf", "both")) {
    ggsave(file.path(plots_dir, "qc_violin_before.pdf"), vln_before, width = 12, height = 4)
    ggsave(file.path(plots_dir, "qc_scatter_before.pdf"), scatter1 + scatter2, width = 10, height = 5)
  }
  if (plot_format %in% c("png", "both")) {
    ggsave(file.path(plots_dir, "qc_violin_before.png"), vln_before, width = 12, height = 4, dpi = 300)
    ggsave(file.path(plots_dir, "qc_scatter_before.png"), scatter1 + scatter2, width = 10, height = 5, dpi = 300)
  }
}

# 5. 执行过滤
message("[R] 执行过滤 (方法: ", cfg$qc_method, ")...")

if (cfg$qc_method != "manual") {
  stop("当前版本仅支持 manual 方法。收到的方法: ", cfg$qc_method)
}

# 手工阈值过滤
pass_qc <- rep(TRUE, ncol(obj))
remove_reasons <- character(ncol(obj))

# 应用阈值
if (cfg$min_features > 0) {
  low_feat <- obj$nFeature_RNA < cfg$min_features
  pass_qc[low_feat] <- FALSE
  remove_reasons[low_feat] <- paste0(remove_reasons[low_feat], ifelse(nchar(remove_reasons[low_feat]) > 0, ";", ""), "low_features")
}

if (cfg$max_features > 0) {
  high_feat <- obj$nFeature_RNA > cfg$max_features
  pass_qc[high_feat] <- FALSE
  remove_reasons[high_feat] <- paste0(remove_reasons[high_feat], ifelse(nchar(remove_reasons[high_feat]) > 0, ";", ""), "high_features")
}

if (cfg$min_counts > 0) {
  low_counts <- obj$nCount_RNA < cfg$min_counts
  pass_qc[low_counts] <- FALSE
  remove_reasons[low_counts] <- paste0(remove_reasons[low_counts], ifelse(nchar(remove_reasons[low_counts]) > 0, ";", ""), "low_counts")
}

if (cfg$max_counts > 0) {
  high_counts <- obj$nCount_RNA > cfg$max_counts
  pass_qc[high_counts] <- FALSE
  remove_reasons[high_counts] <- paste0(remove_reasons[high_counts], ifelse(nchar(remove_reasons[high_counts]) > 0, ";", ""), "high_counts")
}

if (cfg$max_percent_mt > 0) {
  high_mt <- obj[[mt_col_name, drop = TRUE]] > cfg$max_percent_mt
  pass_qc[high_mt] <- FALSE
  remove_reasons[high_mt] <- paste0(remove_reasons[high_mt], ifelse(nchar(remove_reasons[high_mt]) > 0, ";", ""), "high_percent_mt")
}

# 更新 qc_metrics
qc_metrics$pass_qc <- pass_qc
qc_metrics$remove_reason <- remove_reasons

# 6. 生成过滤后对象
message("[R] 生成过滤后对象...")
cells_keep <- colnames(obj)[pass_qc]

if (length(cells_keep) == 0) {
  stop("QC 阈值过严，过滤后无剩余细胞。请调整阈值。")
}

obj_filtered <- subset(obj, cells = cells_keep)
message("[R] 过滤后保留 ", ncol(obj_filtered), " / ", ncol(obj), " 个细胞")

# 7. 保存过滤后对象
output_class <- "Seurat"
if (cfg$output_object == "inherit") {
  output_class <- class(obj)[1]
}

obj_path <- file.path(output_dir, "qc_object.rds")
saveRDS(obj_filtered, obj_path)

# 8. 生成过滤后指标表
qc_kept <- qc_metrics[pass_qc, ]
qc_removed <- qc_metrics[!pass_qc, ]

write.csv(qc_kept, gzfile(file.path(output_dir, "qc_kept_cells.csv.gz")), row.names = FALSE)
write.csv(qc_removed, gzfile(file.path(output_dir, "qc_removed_cells.csv.gz")), row.names = FALSE)

# 9. 生成 QC 后图
if (isTRUE(cfg$save_plots)) {
  message("[R] 生成 QC 后图...")
  
  # Violin plot
  vln_after <- VlnPlot(obj_filtered, features = c("nFeature_RNA", "nCount_RNA", mt_col_name), 
                       pt.size = 0.1, ncol = 3)
  
  # Scatter plots
  scatter1_after <- FeatureScatter(obj_filtered, feature1 = "nCount_RNA", feature2 = mt_col_name)
  scatter2_after <- FeatureScatter(obj_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  # 保存图像
  plot_format <- cfg$plot_format
  if (plot_format %in% c("pdf", "both")) {
    ggsave(file.path(plots_dir, "qc_violin_after.pdf"), vln_after, width = 12, height = 4)
    ggsave(file.path(plots_dir, "qc_scatter_after.pdf"), scatter1_after + scatter2_after, width = 10, height = 5)
  }
  if (plot_format %in% c("png", "both")) {
    ggsave(file.path(plots_dir, "qc_violin_after.png"), vln_after, width = 12, height = 4, dpi = 300)
    ggsave(file.path(plots_dir, "qc_scatter_after.png"), scatter1_after + scatter2_after, width = 10, height = 5, dpi = 300)
  }
}

# 10. 生成统计摘要
remove_counts <- table(unlist(strsplit(remove_reasons[!pass_qc], ";")))
remove_counts <- as.list(remove_counts)

summary <- list(
  tool_name = "qc",
  tool_version = "1.0.0",
  input_object_file = basename(input_path),
  input_object_class = class(obj)[1],
  output_object_file = "qc_object.rds",
  output_object_class = output_class,
  assay = cfg$assay,
  qc_method = cfg$qc_method,
  thresholds = list(
    min_features = cfg$min_features,
    max_features = cfg$max_features,
    min_counts = cfg$min_counts,
    max_counts = cfg$max_counts,
    max_percent_mt = cfg$max_percent_mt,
    mt_pattern = mt_pattern,
    mt_col_name = mt_col_name
  ),
  n_cells_before = ncol(obj),
  n_cells_after = ncol(obj_filtered),
  n_removed = ncol(obj) - ncol(obj_filtered),
  remove_reason_counts = remove_counts,
  plot_files = if (isTRUE(cfg$save_plots)) {
    plot_exts <- if (plot_format == "both") c(".pdf", ".png") else paste0(".", plot_format)
    files <- c()
    for (ext in plot_exts) {
      files <- c(files,
                 paste0("plots/qc_violin_before", ext),
                 paste0("plots/qc_scatter_before", ext),
                 paste0("plots/qc_violin_after", ext),
                 paste0("plots/qc_scatter_after", ext))
    }
    files
  } else character(0)
)

jsonlite::write_json(summary, path = file.path(output_dir, "qc_summary.json"), auto_unbox = TRUE, pretty = TRUE)

message("[R] QC 流程完成")
'''

    @staticmethod
    def _resolve_rscript() -> str:
        """解析正确的 Rscript 路径"""
        import shutil
        import os
        from pathlib import Path

        # 1) 允许手动指定
        for key in ("R_SCRIPT_PATH", "QC_RSCRIPT"):
            value = os.environ.get(key)
            if value and Path(value).exists():
                return value
            
        # 3) 当前 conda 环境
        conda_prefix = os.environ.get("CONDA_PREFIX")
        if conda_prefix:
            candidate = Path(conda_prefix) / "bin" / "Rscript"
            if candidate.exists():
                return str(candidate)

        # 4) 最后才回退到 PATH
        found = shutil.which("Rscript")
        if found:
            return found

        return "Rscript"

    def run(self, script_path: Path, config_path: Path, input_path: Path, output_dir: Path) -> Tuple[int, str, str]:
        """运行 R 脚本"""
        import os
        from pathlib import Path
        
        # 解析正确的 Rscript 路径
        rscript = self._resolve_rscript()
        print(f"[R 环境检查] 使用 Rscript: {rscript}")
        
        # 创建环境变量副本
        env = os.environ.copy()
        rscript_dir = str(Path(rscript).parent)
        env["PATH"] = f"{rscript_dir}:{env.get('PATH', '')}"
        
        # 预检打印：验证 R 环境
        check_cmd = [
            rscript,
            "-e",
            "cat('R.home=', R.home(), '\\n'); print(.libPaths()); library(rlang); library(SeuratObject); library(Seurat); library(ggplot2); cat('R_CHECK_OK\\n')"
        ]
        
        print(f"[R 环境检查] 执行预检命令: {' '.join(check_cmd)}")
        check_proc = subprocess.run(check_cmd, capture_output=True, text=True, env=env)
        if check_proc.returncode != 0:
            print(f"[R 环境检查] 预检失败，返回码: {check_proc.returncode}")
            print(f"[R 环境检查] stdout: {check_proc.stdout}")
            print(f"[R 环境检查] stderr: {check_proc.stderr}")
            # 继续执行主脚本，但记录警告
        else:
            print(f"[R 环境检查] 预检成功")
            print(f"[R 环境检查] R.home: {check_proc.stdout.split('R.home=')[1].splitlines()[0] if 'R.home=' in check_proc.stdout else 'N/A'}")
        
        # 执行主脚本
        cmd = [rscript, str(script_path), str(config_path), str(input_path), str(output_dir)]
        print(f"[R 执行] 执行命令: {' '.join(cmd)}")
        proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
        
        return proc.returncode, proc.stdout, proc.stderr


class ReportWriter:
    @staticmethod
    def write_report(
        report_path: Path,
        config: QCConfig,
        resolved: Dict[str, Any],
        summary: Dict[str, Any],
        stdout_text: str,
        stderr_text: str,
    ) -> None:
        """写入报告文件"""
        
        # 构建阈值描述
        thresholds_desc = []
        if config.min_features > 0:
            thresholds_desc.append(f"- 最小基因数: {config.min_features}")
        if config.max_features > 0:
            thresholds_desc.append(f"- 最大基因数: {config.max_features}")
        if config.min_counts > 0:
            thresholds_desc.append(f"- 最小 UMI 数: {config.min_counts}")
        if config.max_counts > 0:
            thresholds_desc.append(f"- 最大 UMI 数: {config.max_counts}")
        if config.max_percent_mt > 0:
            thresholds_desc.append(f"- 最大线粒体比例: {config.max_percent_mt}%")
        
        thresholds_block = "\n".join(thresholds_desc) if thresholds_desc else "- 无阈值限制"
        
        # 构建删除原因统计
        remove_stats = []
        if "remove_reason_counts" in summary:
            for reason, count in summary["remove_reason_counts"].items():
                remove_stats.append(f"- {reason}: {count}")
        
        remove_stats_block = "\n".join(remove_stats) if remove_stats else "- 无细胞被删除"
        
        report = f"""# 单细胞质控 - 处理报告

## 任务概览
- 项目名称: qc
- 显示名称: 单细胞质控
- 版本: 1.0.0
- 执行时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

## 输入配置
- assay: {config.assay}
- qc_method: {config.qc_method}
- 线粒体基因匹配模式: {config.mt_pattern}
- 线粒体比例列名: {config.mt_col_name}
- 保存图像: {config.save_plots}
- 图像格式: {config.plot_format}
- 输出对象类型: {config.output_object}

## 过滤阈值
{thresholds_block}

## 统计摘要
- 输入对象文件: {summary.get('input_object_file', 'N/A')}
- 输入对象类型: {summary.get('input_object_class', 'N/A')}
- 输出对象文件: {summary.get('output_object_file', 'N/A')}
- 输出对象类型: {summary.get('output_object_class', 'N/A')}
- 过滤前细胞数: {summary.get('n_cells_before', 'N/A')}
- 过滤后细胞数: {summary.get('n_cells_after', 'N/A')}
- 删除细胞数: {summary.get('n_removed', 'N/A')}
- 删除比例: {round(summary.get('n_removed', 0) / max(summary.get('n_cells_before', 1), 1) * 100, 2)}%

## 删除原因统计
{remove_stats_block}

## 输出文件
- qc_object.rds: 过滤后的 Seurat 对象
- qc_metrics_all_cells.csv.gz: 所有细胞的 QC 指标
- qc_kept_cells.csv.gz: 保留细胞的 QC 指标
- qc_removed_cells.csv.gz: 删除细胞的 QC 指标
- qc_summary.json: 统计摘要 JSON
- report.md: 本报告文件

## 图像文件
"""
        
        # 添加图像文件列表
        if config.save_plots and "plot_files" in summary:
            for plot_file in summary["plot_files"]:
                report += f"- {plot_file}\n"
        else:
            report += "- 未生成图像文件\n"
        
        report += f"""
## 运行日志（stdout）
```text
{stdout_text.strip()}
```

## 运行日志（stderr）
```text
{stderr_text.strip()}
```
"""
        report_path.write_text(report, encoding="utf-8")


class runner(Tool):
    """单细胞分析流程中的细胞质量控制步骤"""
    
    DISPLAY_NAME = "单细胞质控"
    NETWORK = False
    CPU = 1
    GPU = 0
    SIF = "singlecellQC:1.0.0"
    
    def call(self, kwargs):
        args = kwargs.get("args", {}) or {}
        config = self._build_config(args)
        
        print("=== 开始执行单细胞 QC 流程 ===")
        print(f"assay: {config.assay}")
        print(f"qc_method: {config.qc_method}")
        print(f"mt_pattern: {config.mt_pattern}")
        print(f"min_features: {config.min_features}")
        print(f"max_features: {config.max_features}")
        print(f"max_percent_mt: {config.max_percent_mt}")
        print(f"save_plots: {config.save_plots}")
        print(f"plot_format: {config.plot_format}")
        
        # 设置临时目录
        cwd = os.getcwd()
        local_tmp_path = os.path.join(cwd, ".tmp")
        os.makedirs(local_tmp_path, exist_ok=True)
        os.environ["TMPDIR"] = local_tmp_path
        
        with tempfile.TemporaryDirectory(prefix="qc_") as tmpdir:
            tmp_path = Path(tmpdir)
            resolver = InputResolver(tmp_path)
            resolved = resolver.resolve(kwargs, config)
            
            print(f"输入文件: {resolved['input_file']}")
            print(f"文件大小: {resolved['file_size']} bytes")
            
            config_path = tmp_path / "config.json"
            output_dir = tmp_path / "outputs"
            config_path.write_text(json.dumps(asdict(config), ensure_ascii=False, indent=2), encoding="utf-8")
            output_dir.mkdir(parents=True, exist_ok=True)
            
            r_builder = RPipelineBuilder(tmp_path)
            script_path = r_builder.write_script()
            returncode, stdout_text, stderr_text = r_builder.run(
                script_path, config_path, Path(resolved["input_path"]), output_dir
            )
            
            print(stdout_text)
            if stderr_text:
                print(stderr_text)
            
            if returncode != 0:
                raise QCError(
                    f"R QC 流程执行失败，返回码: {returncode}。\nstderr:\n{stderr_text}"
                )
            
            # 读取生成的摘要
            summary_path = output_dir / "qc_summary.json"
            if not summary_path.exists():
                raise QCError("QC 流程未生成 qc_summary.json，输出不完整。")
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            
            # 复制输出文件到当前目录
            self._copy_outputs_to_cwd(output_dir)
            
            # 写入报告
            ReportWriter.write_report(
                Path("./report.md"),
                config,
                resolved,
                summary,
                stdout_text,
                stderr_text
            )
        
        print("=== 单细胞 QC 流程执行完成 ===")
        print("已生成文件: report.md, qc_object.rds, qc_*.csv.gz, qc_summary.json")
        if config.save_plots:
            print("已生成图像文件: plots/*")
        
        return json.dumps(
            {
                "status": "success",
                "message": "单细胞 QC 完成",
                "n_cells_before": summary.get("n_cells_before", 0),
                "n_cells_after": summary.get("n_cells_after", 0),
                "n_removed": summary.get("n_removed", 0),
                "output_object": config.output_object,
            },
            ensure_ascii=False,
        )
    
    def outputShow(self, kwargs):
        output = kwargs.get("stdout", "")
        display_text = f"""## ✅ 单细胞质控完成

### 执行摘要
```text
{output}
```

### 输出文件
- `report.md`：详细处理报告
- `qc_object.rds`：过滤后的 Seurat 对象
- `qc_metrics_all_cells.csv.gz`：所有细胞的 QC 指标
- `qc_kept_cells.csv.gz`：保留细胞的 QC 指标
- `qc_removed_cells.csv.gz`：删除细胞的 QC 指标
- `qc_summary.json`：统计摘要 JSON
"""
        
        # 检查是否有图像文件
        import glob
        plot_files = []
        for ext in ["pdf", "png"]:
            plot_files.extend(glob.glob(f"plots/*.{ext}"))
        
        if plot_files:
            display_text += "\n### 图像文件\n"
            for plot_file in sorted(plot_files):
                display_text += f"- `{plot_file}`\n"
        
        file_list = [
            "report.md",
            "qc_object.rds",
            "qc_metrics_all_cells.csv.gz",
            "qc_kept_cells.csv.gz",
            "qc_removed_cells.csv.gz",
            "qc_summary.json",
        ]
        
        # 添加图像文件到文件列表
        file_list.extend(plot_files)
        
        return display_text, file_list
    
    def summary(self, kwargs):
        output = kwargs.get("stdout", "")
        stderr_text = kwargs.get("stderr", "")
        if stderr_text:
            return f"单细胞质控已执行，stdout: {output}；stderr: {stderr_text}"
        return f"单细胞质控已成功完成。{output}"
    
    @staticmethod
    def _build_config(args: Dict[str, Any]) -> QCConfig:
        def _as_bool(value: Any, default: bool) -> bool:
            if isinstance(value, bool):
                return value
            if isinstance(value, str):
                lowered = value.strip().lower()
                if lowered in {"true", "1", "yes", "y"}:
                    return True
                if lowered in {"false", "0", "no", "n"}:
                    return False
            return default
        
        def _as_int(value: Any, default: int) -> int:
            try:
                return int(value)
            except Exception:
                return default
        
        return QCConfig(
            assay=str(args.get("assay", "RNA") or "RNA"),
            qc_method=str(args.get("qc_method", "manual") or "manual"),
            mt_pattern=str(args.get("mt_pattern", "^MT-") or "^MT-"),
            mt_col_name=str(args.get("mt_col_name", "percent.mt") or "percent.mt"),
            min_features=_as_int(args.get("min_features", 200), 200),
            max_features=_as_int(args.get("max_features", 2500), 2500),
            min_counts=_as_int(args.get("min_counts", 0), 0),
            max_counts=_as_int(args.get("max_counts", 0), 0),
            max_percent_mt=_as_int(args.get("max_percent_mt", 5), 5),
            save_plots=_as_bool(args.get("save_plots", True), True),
            plot_format=str(args.get("plot_format", "both") or "both"),
            output_object=str(args.get("output_object", "seurat") or "seurat"),
        )
    
    @staticmethod
    def _copy_outputs_to_cwd(output_dir: Path) -> None:
        """复制输出文件到当前工作目录"""
        expected_files = [
            "qc_object.rds",
            "qc_metrics_all_cells.csv.gz",
            "qc_kept_cells.csv.gz",
            "qc_removed_cells.csv.gz",
            "qc_summary.json",
        ]
        
        # 创建 plots 目录
        plots_dir = Path("./plots")
        plots_dir.mkdir(exist_ok=True)
        
        # 复制主文件
        for filename in expected_files:
            src = output_dir / filename
            if src.exists():
                shutil.copy2(src, Path("./") / filename)
        
        # 复制图像文件
        plots_src_dir = output_dir / "plots"
        if plots_src_dir.exists():
            for plot_file in plots_src_dir.glob("*"):
                shutil.copy2(plot_file, plots_dir / plot_file.name)
