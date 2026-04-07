parent_path <- ifelse(Sys.info()["sysname"] == "Windows", "G:", "/nfs/user/Users/wl2907");
project_ID <- "PROJ0111";
data_label <- "KIF1A";

load_data <- TRUE;
force_run <- TRUE;

# === Load Saved .RData ===
if (load_data) {
    filename <- dir(paste(parent_path, "/Projects/", project_ID, sep=""), full.names = TRUE);
    filename <- filename[grepl(".RData$", filename) & grepl(data_label, filename)];
    if (length(filename) > 0) {
        filename <- sort(filename, decreasing = TRUE)[1];
        load(filename);
        message(paste("Data: \"", filename, "\" loaded.", sep=""));
    } else {
        message(paste("Data: \"", filename, "\" not found.", sep=""));
    }
    parent_path <- ifelse(Sys.info()["sysname"] == "Windows", "G:", "/nfs/user/Users/wl2907");
}

project_path <- paste(parent_path, "/Projects/", project_ID, sep="");
resource_path <- paste(parent_path, "/Resources", sep="");
tool_path <- paste(parent_path, "/Tools", sep="");

# Data I/O
library(openxlsx);       # Read/write Excel files
library(rentrez);        # Fetch NCBI sequence records
library(httr);           # HTTP requests (REST API calls)
library(jsonlite);       # JSON parsing
# Data Wrangling
library(dplyr);          # Data manipulation
# Bioinformatics
library(Biostrings);     # DNA/protein sequence manipulation
library(seqinr);         # Amino acid code conversion
# Statistics & Machine Learning
library(lcmm);           # Group-based trajectory modeling
library(glmnet);         # Regularized regression (LASSO/Ridge)
library(randomForest);   # Random forest models
library(caret);          # ML training framework
library(pROC);           # ROC curve analysis
library(factoextra);     # PCA visualization (fviz_eig)
# Visualization: General
library(ggplot2);        # Core plotting
library(ggrepel);        # Non-overlapping text labels
library(ggpmisc);        # Stat annotations on ggplot
library(ggpubr);         # Publication-ready plots (ggarrange, stat_compare_means)
library(ggcorrplot);     # Correlation matrix plots
library(ggdendro);       # Dendrogram visualization
library(plotly);         # Interactive plots
library(cowplot);        # Plot composition and theming
library(patchwork);      # Combine ggplot panels
# Visualization: Heatmap
library(pheatmap);       # Static heatmaps
library(RColorBrewer);   # Color palettes for heatmaps
library(ComplexHeatmap); # Advanced heatmaps with annotations
library(circlize);       # Circular plots and ColorRamp2 for heatmaps
# Visualization: Layout & Rendering
library(grid);           # Low-level graphics grid system
library(gridExtra);      # Arrange multiple grid-based plots
library(htmlwidgets);    # Export interactive widgets to HTML
# Computing
library(parallel);       # Parallel computing (mclapply, makeCluster)

setwd(project_path);
source(paste(tool_path, "/General Drawing Functions.R", sep=""));
source(paste(tool_path, "/Genetics/De novo Functions.R", sep=""));
source(paste(tool_path, "/Clinical/Clinical Data Analysis Functions.R", sep=""));

rawdata_path <- paste(project_path, "/1_Raw_data", sep="");

# === Function Start ===
get_codon_info <- function(transcript_id, pos, alt_base) {
    # Fetch GenBank record to extract CDS position
    gb <- entrez_fetch(db = "nuccore",
                        id = transcript_id,
                        rettype = "gb",
                        retmode = "text");
    
    # Parse CDS start position
    cds_match <- regmatches(gb, regexpr("CDS\\s+\\d+\\.\\.\\d+", gb));
    cds_start <- as.integer(regmatches(cds_match, regexpr("\\d+", cds_match)));
    
    # Fetch FASTA sequence
    fasta <- entrez_fetch(db = "nuccore",
                            id = transcript_id,
                            rettype = "fasta",
                            retmode = "text");
    
    # Parse sequence string
    seq_lines <- strsplit(fasta, "\n")[[1]];
    seq_str <- paste(seq_lines[-1], collapse = "");
    rna_seq <- DNAString(seq_str);
    
    # Calculate absolute position in transcript
    abs_pos <- cds_start + pos - 1;
    
    # Determine codon information
    codon_num   <- ceiling(pos / 3);              # Amino acid number
    codon_pos   <- ((pos - 1) %% 3) + 1;         # Position within codon (1/2/3)
    codon_start <- cds_start + (codon_num - 1) * 3;  # Codon start in transcript
    
    # Extract reference base and codon
    ref_base <- as.character(rna_seq[abs_pos]);
    codon    <- rna_seq[codon_start:(codon_start + 2)];
    aa_ref   <- as.character(translate(codon));
    
    # Construct mutant codon and translate
    codon_mut <- codon;
    codon_mut[codon_pos] <- alt_base;
    aa_mut <- as.character(translate(codon_mut));
    
    # Determine mutation consequence
    if (ref_base == alt_base) {
        consequence <- "no change";
    } else if (aa_ref == aa_mut) {
        consequence <- "synonymous";
    } else if (aa_mut == "*") {
        consequence <- "nonsense (stop gained)";
    } else if (aa_ref == "*") {
        consequence <- "stop lost";
    } else {
        consequence <- "missense";
    };
    
    # Print results
    cat("─────────────────────────────────\n");
    cat("Transcript:         ", transcript_id, "\n");
    cat("c.position:         ", pos, "\n");
    cat("CDS start:          ", cds_start, "\n");
    cat("─────────────────────────────────\n");
    cat("Reference base:     ", ref_base, "\n");
    cat("Codon number:       ", codon_num, "\n");
    cat("Position in codon:  ", codon_pos, "\n");
    cat("Reference codon:    ", as.character(codon), "->", aa_ref, "\n");
    cat("─────────────────────────────────\n");
    cat("Alt base:           ", alt_base, "\n");
    cat("Mutant codon:       ", as.character(codon_mut), "->", aa_mut, "\n");
    cat("Consequence:        ", consequence, "\n");
    cat("HGVS (c.):         ",
        paste0("c.", pos, ref_base, ">", alt_base), "\n");
    cat("HGVS (p.):         ",
        paste0("p.", aa_ref, codon_num, aa_mut), "\n");
    cat("─────────────────────────────────\n");
    
    # Return results as a list
    invisible(list(
        transcript         = transcript_id,
        pos                = pos,
        cds_start          = cds_start,
        abs_pos            = abs_pos,
        ref_base           = ref_base,
        alt_base           = alt_base,
        codon_num          = codon_num,
        codon_pos_in_codon = codon_pos,
        ref_codon          = as.character(codon),
        mut_codon          = as.character(codon_mut),
        aa_ref             = aa_ref,
        aa_mut             = aa_mut,
        consequence        = consequence
    ));
}

PCA_analysis <- function(clinical_data, out_path = NULL, k_max = 8, best_k = NULL) {
    if (is.null(out_path)) out_path <- getwd();
    if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE);

    # ── 1. Convert to numeric and remove columns with NA / NaN ───────────────
    for (i in seq_len(ncol(clinical_data))) {
        clinical_data[, i] <- as.numeric(clinical_data[, i]);
    }
    clinical_data <- clinical_data[, !apply(clinical_data, 2, function(x) any(is.na(x))), drop = FALSE];

    data_normalized <- scale(clinical_data);
    data_normalized <- data_normalized[, !apply(data_normalized, 2, function(x) any(is.nan(x))), drop = FALSE];

    # ── 2. PCA ────────────────────────────────────────────────────────────────
    data.pca <- prcomp(data_normalized);
    print(summary(data.pca));

    var_exp <- summary(data.pca)$importance["Proportion of Variance", ];
    pc1_pct <- round(var_exp["PC1"] * 100, 1);
    pc2_pct <- round(var_exp["PC2"] * 100, 1);

    scores       <- as.data.frame(data.pca$x[, 1:2]);
    scores$label <- rownames(clinical_data);

    # ── 3. White background theme ─────────────────────────────────────────────
    theme_white <- ggplot2::theme_classic(base_size = 12) +
        ggplot2::theme(
            plot.background   = ggplot2::element_rect(fill = "white", color = NA),
            panel.background  = ggplot2::element_rect(fill = "white", color = NA),
            legend.background = ggplot2::element_rect(fill = "white", color = NA),
            plot.title        = ggplot2::element_text(size = 12, face = "plain")
        );

    theme_white_patch <- ggplot2::theme(
        plot.background   = ggplot2::element_rect(fill = "white", color = NA),
        panel.background  = ggplot2::element_rect(fill = "white", color = NA),
        legend.background = ggplot2::element_rect(fill = "white", color = NA)
    );

    # ── 4. Legacy plots with white background ─────────────────────────────────
    p <- ggcorrplot::ggcorrplot(cor(data_normalized)) + theme_white_patch;
    ggplot2::ggsave(file.path(out_path, "Correlation matrix.png"), p, width = 12, height = 10, bg = "white");

    p <- factoextra::fviz_eig(data.pca, addlabels = TRUE, ncp = 20) + theme_white_patch;
    ggplot2::ggsave(file.path(out_path, "Scree Plot.png"), p, width = 8, height = 8, bg = "white");

    p <- factoextra::fviz_pca_var(data.pca, col.var = "black") + theme_white_patch;
    ggplot2::ggsave(file.path(out_path, "Biplot of the attributes.png"), p, width = 12, height = 8, bg = "white");

    for (i in 1:3) {
        p <- factoextra::fviz_cos2(data.pca, choice = "var", axes = i) +
            theme_white_patch +
            ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 1, "cm"));
        ggplot2::ggsave(
            file.path(out_path, paste0("Contribution of each variable on PC", i, ".png")),
            p, width = 8, height = 8, bg = "white"
        );
    }

    p <- factoextra::fviz_pca_var(data.pca, col.var = "cos2",
                                  gradient.cols = c("black", "orange", "green"), repel = TRUE) +
        theme_white_patch;
    ggplot2::ggsave(file.path(out_path, "Biplot combined with cos2.png"), p, width = 12, height = 8, bg = "white");

    # ── 5. Variance explained bar chart ───────────────────────────────────────
    n_pc   <- min(k_max, length(var_exp));
    df_var <- data.frame(
        PC       = factor(paste0("PC", seq_len(n_pc)), levels = paste0("PC", seq_len(n_pc))),
        Variance = as.numeric(var_exp[seq_len(n_pc)]) * 100
    );

    p <- ggplot2::ggplot(df_var, ggplot2::aes(x = PC, y = Variance)) +
        ggplot2::geom_col(fill = "#4472C4", width = 0.6) +
        ggplot2::geom_text(ggplot2::aes(label = paste0(round(Variance, 1), "%")),
                           vjust = -0.4, size = 3.5) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
        ggplot2::labs(title = "Variance Explained by PCs", x = "PC", y = "Variance (%)") +
        theme_white;
    ggplot2::ggsave(file.path(out_path, "Variance Explained by PCs.png"), p,
                    width = 6, height = 5, bg = "white");

    # ── 6. k-means clustering + silhouette + elbow + gap statistic ───────────
    scores_mat <- as.matrix(scores[, c("PC1", "PC2")]);

    sil_avg <- sapply(2:k_max, function(k) {
        km  <- kmeans(scores_mat, centers = k, nstart = 25, iter.max = 100);
        sil <- cluster::silhouette(km$cluster, dist(scores_mat));
        mean(sil[, 3]);
    });

    wss_vals <- sapply(2:k_max, function(k) {
        km <- kmeans(scores_mat, centers = k, nstart = 25, iter.max = 100);
        km$tot.withinss;
    });

    # Gap statistic (used when best_k is NULL)
    if (is.null(best_k)) {
        message("Computing Gap statistic ...");
        gap_result <- cluster::clusGap(
            scores_mat,
            FUN     = kmeans,
            K.max   = k_max,
            B       = 100,
            verbose = FALSE,
            nstart  = 25,
            iter.max = 100
        );
        best_k <- cluster::maxSE(
            gap_result$Tab[, "gap"],
            gap_result$Tab[, "SE.sim"],
            method = "Tibs2001SEmax"
        );
        message("Optimal k by Gap statistic: ", best_k);

        # Gap statistic plot
        gap_df <- data.frame(
            k   = seq_len(k_max),
            gap = gap_result$Tab[, "gap"],
            se  = gap_result$Tab[, "SE.sim"]
        );

        p <- ggplot2::ggplot(gap_df, ggplot2::aes(x = k, y = gap)) +
            ggplot2::geom_line(color = "#4472C4", linewidth = 0.8) +
            ggplot2::geom_point(size = 3, color = "#4472C4") +
            ggplot2::geom_errorbar(ggplot2::aes(ymin = gap - se, ymax = gap + se),
                                   width = 0.2, color = "#4472C4") +
            ggplot2::geom_vline(xintercept = best_k,
                                color = "red", linetype = "dashed", linewidth = 0.7) +
            ggplot2::scale_x_continuous(breaks = seq_len(k_max)) +
            ggplot2::labs(title = "Gap statistic",
                          x = "Number of clusters k", y = "Gap statistic") +
            theme_white;
        ggplot2::ggsave(file.path(out_path, "Gap statistic.png"), p,
                        width = 6, height = 5, bg = "white");
    };

    message("Using k = ", best_k, " for clustering");

    # Silhouette plot
    df_sil <- data.frame(k = 2:k_max, silhouette = sil_avg);

    p <- ggplot2::ggplot(df_sil, ggplot2::aes(x = k, y = silhouette)) +
        ggplot2::geom_line(color = "#4472C4", linewidth = 0.8) +
        ggplot2::geom_point(size = 3, color = "#4472C4") +
        ggplot2::geom_point(data = df_sil[which.max(df_sil$silhouette), ],
                            size = 5, shape = 21,
                            fill = "#4472C4", color = "white", stroke = 1.5) +
        ggplot2::geom_vline(xintercept = best_k,
                            color = "red", linetype = "dashed", linewidth = 0.7) +
        ggplot2::scale_x_continuous(breaks = 2:k_max) +
        ggplot2::labs(title = "Silhouette plot",
                      x = "Number of clusters k", y = "Average silhouette width") +
        theme_white;
    ggplot2::ggsave(file.path(out_path, "Silhouette plot.png"), p,
                    width = 6, height = 5, bg = "white");

    # Elbow plot
    df_wss <- data.frame(k = 2:k_max, wss = wss_vals);

    p <- ggplot2::ggplot(df_wss, ggplot2::aes(x = k, y = wss)) +
        ggplot2::geom_line(color = "#4472C4", linewidth = 0.8) +
        ggplot2::geom_point(size = 3, color = "#4472C4") +
        ggplot2::geom_vline(xintercept = best_k,
                            color = "red", linetype = "dashed", linewidth = 0.7) +
        ggplot2::scale_x_continuous(breaks = 2:k_max) +
        ggplot2::labs(title = "Elbow plot",
                      x = "Number of clusters k", y = "Total within-cluster SS") +
        theme_white;
    ggplot2::ggsave(file.path(out_path, "Elbow plot.png"), p,
                    width = 6, height = 5, bg = "white");

    # ── 7. Hierarchical clustering dendrogram ─────────────────────────────────
    hc         <- hclust(dist(scores_mat), method = "ward.D2");
    cut_height <- mean(c(hc$height[length(hc$height) - best_k + 1],
                         hc$height[length(hc$height) - best_k + 2]));

    p <- ggdendro::ggdendrogram(hc, rotate = FALSE, labels = FALSE) +
        ggplot2::geom_hline(yintercept = cut_height,
                            color = "red", linetype = "dashed", linewidth = 0.7) +
        ggplot2::labs(title = "Dendrogram (ward.D2)") +
        ggplot2::theme(
            plot.background   = ggplot2::element_rect(fill = "white", color = NA),
            panel.background  = ggplot2::element_rect(fill = "white", color = NA),
            plot.title        = ggplot2::element_text(size = 12, face = "plain"),
            axis.text.x       = ggplot2::element_blank()
        );
    ggplot2::ggsave(file.path(out_path, "Dendrogram.png"), p,
                    width = 6, height = 5, bg = "white");

    # ── 8. PCA scatter plot (k-means coloring) + row name labels ─────────────
    km_best        <- kmeans(scores_mat, centers = best_k, nstart = 25, iter.max = 100);
    scores$Cluster <- factor(km_best$cluster);

    cluster_colors <- setNames(
        scales::hue_pal()(best_k),
        as.character(seq_len(best_k))
    );

    p <- ggplot2::ggplot(scores, ggplot2::aes(x = PC1, y = PC2, color = Cluster)) +
        ggplot2::geom_point(size = 3, alpha = 0.85) +
        ggrepel::geom_text_repel(
            ggplot2::aes(label = label),
            size          = 2.8,
            max.overlaps  = Inf,
            segment.size  = 0.3,
            segment.alpha = 0.5,
            show.legend   = FALSE
        ) +
        ggplot2::scale_color_manual(values = cluster_colors, name = "Cluster") +
        ggplot2::labs(
            x = paste0("PC 1 (", pc1_pct, "%)"),
            y = paste0("PC 2 (", pc2_pct, "%)")
        ) +
        theme_white +
        ggplot2::theme(legend.position = "right");
    ggplot2::ggsave(file.path(out_path, "PC1 vs PC2 clustering.png"),
                    p, width = 10, height = 8, dpi = 300, bg = "white");

    # ── 9. Export PCA_result.xlsx ─────────────────────────────────────────────
    pc_scores <- as.data.frame(data.pca$x);
    colnames(pc_scores) <- paste0("PC", seq_len(ncol(pc_scores)));

    pca_cluster <- km_best$cluster[match(rownames(clinical_data), rownames(scores))];

    pca_result <- cbind(clinical_data, pc_scores, PCA_cluster = pca_cluster);
    pca_result <- data.frame(ID = rownames(pca_result), pca_result, check.names = FALSE);

    openxlsx::write.xlsx(pca_result, file.path(out_path, "PCA_result.xlsx"),
                         rowNames = FALSE);

    return(list(pca = data.pca, result = pca_result));
};

calc_dim_reduction <- function(data, preprocessing=TRUE, add_noise=TRUE, method="PCA") {
    if (preprocessing == TRUE) {
        # Convert data to numeric
        for (i in c(1:ncol(data))) {
            data[, i] <- as.numeric(data[, i]);
        }
        data <- data[, !apply(data, 2, function(x) { return(any(is.na(x))); })];
        # Normalizing the data
        data.scaled <- scale(data);
        data.scaled <- data.scaled[, !apply(data.scaled, 2, function(x) { return(any(is.nan(x))); })];
    } else {
        data.scaled <- data;
    }
    # Add a micro noise matrix to the data
    if (add_noise == TRUE) {
        data.scaled <- data.scaled + matrix(rnorm(prod(dim(data.scaled)), mean=0, sd=1e-6), nrow=nrow(data.scaled));
    }
    # Perform unsupervised clustering
    if (tolower(method) == "pca") {
        data.processed <- prcomp(data.scaled);
    } else if (tolower(method) == "tsne") {
        cores <- as.integer(detectCores() / 2);
        data.processed <- Rtsne(data.scaled, pca=FALSE, perplexity=30, theta=0.0, num_threads=cores);
        rownames(data.processed$Y) <- rownames(data);
        colnames(data.processed$Y) <- c("tSNE_1", "tSNE_2");
    } else if (tolower(method) == "umap") {
        data.processed <- umap(data.scaled);
        colnames(data.processed$layout) <- c("UMAP_1", "UMAP_2");
    } else {
        stop("Parameter method must be \"PCA\", \"t-SNE\", or \"UMAP\".");
    }
    # Return the results
    return(data.processed);
}

get_rr_p_label <- function(rr, p, n = 1, rr.digits = 2, p.digits = 3) {
    if (rr < 0.01) {
        rr_label <- "`<`~0.01";
    } else {
        rr_label <- paste("`=`~\"", formatC(rr, format = "f", digits = rr.digits), "\"", sep="");
    }
    if (n == 1) {
        if (p < 0.001) {
            p_label <- "`<`~0.001";
        } else {
            p_label <- paste("`=`~\"", formatC(p, format = "f", digits = p.digits), "\"", sep="");
        }
        rr_p_label <- paste("italic(R)^2~", rr_label, "*`,`~italic(P)~", p_label, sep="");
    } else if (n > 1) {
        p_adj <- p * n;
        if (p_adj < 0.001) {
            p_label <- "`<`~0.001";
        } else if (p_adj >= 1.0) {
            p_label <- "`>`~0.999";
        } else {
            p_label <- paste("`=`~\"", formatC(p_adj, format = "f", digits = p.digits), "\"", sep="");
        }
        rr_p_label <- paste("italic(R)^2~", rr_label, "*`,`~italic(P)[adj]~", p_label, sep="");
    }
    return(rr_p_label);
}

VABS_trend.line_plot <- function(data, anno, VABS_name, var_name = NULL, color_list = NULL, 
    x_lim = NULL, y_lim = NULL, draw_legend = TRUE, height = 4, width = 6, filename = NULL) {
    
    # Define instance values (up to 5 time points)
    instance_vals <- sort(unique(data$Instance));
    
    # Get age and VABS score for each record/instance combination
    age_col <- "Age_VABS";
    VABS_col <- VABS_name;
    
    # Compute x and y limits from long-format data
    if (is.null(x_lim)) {
        age_list <- data[[age_col]];
        age_list <- age_list[!is.na(age_list)];
        x_lim <- range(age_list);
    }
    if (is.null(y_lim)) {
        VABS_list <- data[["VABS_ABC"]];
        VABS_list <- VABS_list[!is.na(VABS_list)];
        y_lim <- range(VABS_list);
    }
    if (is.null(filename)) {
        filename <- paste(VABS_name, "_trend.line_plot.png", sep = "");
    }
    
    # Build color_data from unique Record_IDs
    if (!is.null(var_name)) {
        # Use one row per patient for the coloring variable
        color_data <- unique(data[, c("Record_ID", var_name)]);
        color_data$Color <- NA;
        var_type <- get_var_type(anno, var_name);
        if (var_type == 1) {
            var_list <- quantile(color_data[, var_name], na.rm = TRUE);
            if (is.null(color_list) | length(color_list) < length(var_list)) {
                color_list <- rainbow(length(var_list));
            }
            for (i in c(1:nrow(color_data))) {
                if (!is.na(color_data[i, var_name])) {
                    color_data$Color[i] <- gradient_color(color_data[i, var_name], var_list, color_list);
                }
            }
        } else if (var_type == 2) {
            var_anno <- get_var_annotation(anno, var_name);
            var_list <- var_anno$Detail;
            if (is.null(color_list) | length(color_list) < length(var_list)) {
                color_list <- rainbow(length(var_list));
            }
            color_list <- color_list[1:length(var_list)];
            names(color_list) <- var_list;
            for (i in c(1:nrow(color_data))) {
                if (!is.na(color_data[i, var_name])) {
                    temp_value <- var_anno[var_anno$Item == color_data[i, var_name], ]$Detail;
                    color_data$Color[i] <- color_list[temp_value];
                }
            }
        } else {
            var_freq <- as.data.frame(table(data[, var_name], dnn = var_name));
            var_freq <- var_freq[order(var_freq$Freq, decreasing = TRUE), ];
            var_freq[, var_name] <- as.character(var_freq[, var_name]);
            var_freq <- var_freq[!is.na(var_freq[, var_name]) & var_freq[, var_name] != ".", ];
            count <- min(10, nrow(var_freq));
            if (is.null(color_list) | length(color_list) < count) {
                color_list <- rainbow(count);
            }
            color_list <- color_list[1:count];
            var_list <- var_freq[1:count, var_name];
            names(color_list) <- var_list;
            for (i in c(1:nrow(color_data))) {
                if (!is.na(color_data[i, var_name]) & color_data[i, var_name] %in% names(color_list)) {
                    color_data$Color[i] <- color_list[as.character(color_data[i, var_name])];
                }
            }
        }
    } else {
        color_data <- data.frame(Record_ID = unique(data$Record_ID), Color = NA, stringsAsFactors = FALSE);
    }
    
    # Helper: extract ordered x/y vectors for a single patient from long-format data
    get_xy <- function(ID) {
        temp <- data[data$Record_ID == ID, ];
        temp <- temp[order(temp$Instance), ];
        x <- temp[[age_col]];
        y <- temp[[VABS_col]];
        # Keep only paired non-NA observations
        valid <- !is.na(x) & !is.na(y);
        list(x = x[valid], y = y[valid]);
    }
    
    # Plotting
    png(filename, width = width, height = height, units = "in", res = 600);
    if (draw_legend) {
        layout(matrix(c(rep(1, 5), 2), nrow = 1));
    }
    par(mar = c(3.5, 3.5, 1, 1), mgp = c(2.0, 0.6, 0));
    plot(NULL, NULL, xlim = x_lim, ylim = y_lim, xlab = "Age at evaluation (years)", 
         ylab = paste(gsub("_", " ", VABS_name), "Score"), las = 1);
    
    # Draw lines and points for patients without assigned color (gray background)
    unique_IDs <- color_data[is.na(color_data$Color), ]$Record_ID;
    if (length(unique_IDs) > 0) {
        for (ID in unique_IDs) {
            temp_data <- data[data$Record_ID == ID, ];
            sex_vals <- temp_data[!is.na(temp_data$Sex), ]$Sex;
            sex <- if (length(sex_vals) > 0) sex_vals[1] else NA;
            shape <- ifelse(!is.na(sex) & sex == 1, 16, 17);
            xy <- get_xy(ID);
            if (length(xy$x) > 1) lines(xy$x, xy$y, cex = 0.8, col = "#C5C5C5");
            if (length(xy$x) > 0) points(xy$x, xy$y, pch = shape, col = "#C5C5C5");
        }
    }
    
    # Draw lines and points for patients with assigned colors
    unique_IDs <- color_data[!is.na(color_data$Color), ]$Record_ID;
    if (length(unique_IDs) > 0) {
        for (ID in unique_IDs) {
            temp_data <- data[data$Record_ID == ID, ];
            sex_vals <- temp_data[!is.na(temp_data$Sex), ]$Sex;
            sex <- if (length(sex_vals) > 0) sex_vals[1] else NA;
            shape <- ifelse(!is.na(sex) & sex == 1, 16, 17);
            color <- color_data[color_data$Record_ID == ID, ]$Color;
            xy <- get_xy(ID);
            if (length(xy$x) > 1) lines(xy$x, xy$y, lwd = 1.5, col = color);
            if (length(xy$x) > 0) points(xy$x, xy$y, pch = shape, col = color);
        }
    }
    
    # Draw legend panel
    if (draw_legend) {
        blank_figure(c(0, 10), c(-3, 15), FALSE);
        # Sex legend
        text(0, 14, "Sex", pos = 4);
        points(1, 13, pch = 16, cex = 1.5);
        text(1.2, 13, "Male", pos = 4);
        points(1, 12, pch = 17, cex = 1.5);
        text(1.2, 12, "Female", pos = 4);
        # Variable color legend
        if (!is.null(var_name)) {
            text(0, 10.5, var_name, pos = 4);
            var_type <- get_var_type(anno, var_name);
            if (var_type == 1) {
                temp_value <- c(seq(var_list[5], var_list[4], length.out = 25), 
                                seq(var_list[4], var_list[3], length.out = 25), 
                                seq(var_list[3], var_list[2], length.out = 25), 
                                seq(var_list[2], var_list[1], length.out = 25));
                temp_color <- sapply(temp_value, function(x) { gradient_color(x, var_list, color_list); });
                y_pos <- seq(9.5, 5.5, length.out = 101);
                rect(xleft = 0.8, ybottom = y_pos[2:101], xright = 2, ytop = y_pos[1:100], 
                    col = temp_color, border = NA);
                for (i in c(1:5)) {
                    text(2, 10.5 - i, sprintf("%0.1f", var_list[6 - i]), pos = 4);
                }
            } else {
                for (i in c(1:length(color_list))) {
                    points(1, 10.5 - i, pch = 15, cex = 1.5, col = color_list[i]);
                    text(1.2, 10.5 - i, names(color_list)[i], pos = 4);
                }
            }
        }
    }
    temp <- dev.off();
}

VABS_trend.scatter_plot <- function(KIF1A_merge, variant_freq, min_freq = NULL, top_n = NULL,
                                    colors_list = NULL, show_legend = TRUE,
                                    height = 5, width = 8, filename = NULL) {
    # Determine variant list: min_freq takes priority over top_n
    if (!is.null(min_freq)) {
        variant_list <- variant_freq[variant_freq$Freq >= min_freq, ]$Variant_name;
    } else if (!is.null(top_n)) {
        variant_freq_sorted <- variant_freq[order(variant_freq$Freq, decreasing = TRUE), ];
        variant_list <- variant_freq_sorted$Variant_name[1:min(top_n, nrow(variant_freq_sorted))];
    } else {
        stop("Either min_freq or top_n must be specified.");
    }

    # Filter data: remove NAs, keep valid Instance, keep selected variants
    KIF1A_temp <- KIF1A_merge[
        !is.na(KIF1A_merge$Variant_name) &
        !is.na(KIF1A_merge$Age_VABS) &
        !is.na(KIF1A_merge$VABS_ABC) &
        KIF1A_merge$Instance > 0 &
        KIF1A_merge$Variant_name %in% variant_list, ];

    # Set factor levels to control color assignment order
    KIF1A_temp$Variant_name <- factor(KIF1A_temp$Variant_name, levels = variant_list);

    # Default color palette if not provided
    if (is.null(colors_list)) {
        colors_list <- rainbow(length(variant_list));
    }

    # Build plot
    p <- ggplot(KIF1A_temp, aes(x = Age_VABS, y = VABS_ABC)) +
        geom_point(aes(color = Variant_name), size = 2, alpha = 0.5) +
        geom_smooth(aes(color = Variant_name, group = Variant_name), method = "loess", se = FALSE,
                    linewidth = 1, span = 10) +
        labs(
            x = "Age at evaluation (years)",
            y = "VABS ABC score",
            color = "Variant"
        ) +
        scale_color_manual(values = colors_list) +
        theme_classic() + 
        theme(
            legend.position = ifelse(show_legend, "right", "none"),
            axis.title = element_text(color = "black"),   # Axis titles
            axis.text  = element_text(color = "black"),   # Axis tick labels
            axis.ticks = element_line(color = "black")    # Axis ticks
        );

    # Return plot if no filename specified; otherwise save to file
    if (is.null(filename)) {
        return(p);
    } else {
        ggsave(filename, plot = p, height = height, width = width, units = "in", dpi = 600);
    }
}

VABS_trend.cluster <- function(
        data,
        id_col           = "Record_ID",
        instance_col     = "Instance",
        age_col          = "Age_VABS",
        vabs_col         = "VABS_ABC",
        spline_df        = 3,
        k                = NULL,
        mean_vabs_weight = 2,
        slope_weight     = 1,
        seed             = 42,
        color_list       = NULL,
        filepath         = NULL,
        plot_width       = 8,
        plot_height      = 6,
        plot_dpi         = 300
    ) {

    # ── 0. Helper: save a ggplot to filepath if specified ─────────────────────
    .save_plot <- function(plot_obj, filename) {
        if (!is.null(filepath)) {
            if (!dir.exists(filepath)) {
                dir.create(filepath, recursive = TRUE);
                message("Created directory: ", filepath);
            };
            full_path <- file.path(filepath, filename);
            ggplot2::ggsave(
                filename = full_path,
                plot     = plot_obj,
                width    = plot_width,
                height   = plot_height,
                dpi      = plot_dpi
            );
            message("Saved: ", full_path);
        };
    };

    # ── 1. Rename and clean columns ───────────────────────────────────────────
    df <- data %>%
        dplyr::rename(
            Record_ID = !!rlang::sym(id_col),
            Instance  = !!rlang::sym(instance_col),
            Age_VABS  = !!rlang::sym(age_col),
            VABS_ABC  = !!rlang::sym(vabs_col)
        ) %>%
        dplyr::filter(Instance > 0, !is.na(VABS_ABC), !is.na(Age_VABS)) %>%
        dplyr::arrange(Record_ID, Instance);

    # ── 2. Fit LME model with natural spline (residuals for reference only) ───
    lme_formula <- stats::as.formula(
        paste0("VABS_ABC ~ splines::ns(Age_VABS, df = ", spline_df, ") + (1 | Record_ID)")
    );

    lme_model   <- lme4::lmer(lme_formula, data = df, REML = TRUE);
    df$residual <- stats::residuals(lme_model);

    # ── 3. Extract per-patient clustering features ────────────────────────────
    patient_features <- df %>%
        dplyr::group_by(Record_ID) %>%
        dplyr::summarise(
            n_instance    = dplyr::n(),
            baseline_vabs = dplyr::first(VABS_ABC),
            mean_vabs     = mean(VABS_ABC),
            age_start     = dplyr::first(Age_VABS),
            age_mean      = mean(Age_VABS),
            slope_vabs    = ifelse(
                dplyr::n() >= 2,
                stats::coef(stats::lm(VABS_ABC ~ Age_VABS))[2],
                NA_real_
            ),
            .groups = "drop"
        );

    single_pts   <- patient_features %>% dplyr::filter(n_instance == 1);
    cluster_data <- patient_features %>% dplyr::filter(n_instance >= 2);

    n_clustered <- nrow(cluster_data);
    n_single    <- nrow(single_pts);

    message("Patients eligible for clustering (>= 2 instances): ", n_clustered);
    message("Patients with single instance (excluded from clustering): ", n_single);

    # ── 4. Scale features and apply weights ───────────────────────────────────
    mean_vabs_scaled <- scale(cluster_data$mean_vabs)  * mean_vabs_weight;
    slope_scaled     <- scale(cluster_data$slope_vabs) * slope_weight;

    features_scaled  <- cbind(mean_vabs_scaled, slope_scaled);
    colnames(features_scaled) <- c("mean_vabs", "slope_vabs");

    # ── 5. Diagnostic plots ───────────────────────────────────────────────────
    set.seed(seed);
    plot_elbow <- factoextra::fviz_nbclust(
        features_scaled, kmeans, method = "wss", k.max = min(8, n_clustered - 1)
    ) + ggplot2::labs(title = "Elbow Method – Optimal Number of Clusters");

    plot_silhouette <- factoextra::fviz_nbclust(
        features_scaled, kmeans, method = "silhouette", k.max = min(8, n_clustered - 1)
    ) + ggplot2::labs(title = "Silhouette Method – Optimal Number of Clusters");

    .save_plot(plot_elbow,      "diagnostic_elbow.png");
    .save_plot(plot_silhouette, "diagnostic_silhouette.png");

    # ── 6. Early return if k not specified ────────────────────────────────────
    if (is.null(k)) {
        message("k not specified. Diagnostic plots saved. Set k to proceed with clustering.");
        return(list(
            model               = lme_model,
            data_with_residuals = df,
            patient_features    = patient_features,
            n_clustered         = n_clustered,
            n_single            = n_single
        ));
    };

    # ── 7. K-means clustering ─────────────────────────────────────────────────
    set.seed(seed);
    km_result                <- kmeans(features_scaled, centers = k, nstart = 25);
    cluster_data$cluster_raw <- km_result$cluster;

    # Re-label by ascending mean_vabs: Cluster 1 = lowest (most severe)
    mean_vabs_per_cluster <- tapply(
        cluster_data$mean_vabs,
        cluster_data$cluster_raw,
        mean
    );
    vabs_rank            <- rank(mean_vabs_per_cluster);
    recode_map           <- setNames(as.character(vabs_rank), names(vabs_rank));
    cluster_data$cluster <- factor(
        recode_map[as.character(cluster_data$cluster_raw)],
        levels = as.character(seq_len(k))
    );

    # ── 8. Merge cluster labels back to full data ─────────────────────────────
    data_plot <- df %>%
        dplyr::left_join(
            cluster_data %>% dplyr::select(Record_ID, cluster),
            by = "Record_ID"
        ) %>%
        dplyr::mutate(
            cluster = dplyr::if_else(is.na(cluster), "Single", as.character(cluster)),
            cluster = factor(cluster, levels = c("Single", as.character(seq_len(k))))
        );

    # ── 9. Define color palette ───────────────────────────────────────────────
    if (is.null(color_list)) {
        cluster_colors <- stats::setNames(
            c(grDevices::rainbow(k), "grey70"),
            c(as.character(seq_len(k)), "Single")
        );
    } else {
        if (length(color_list) < k) {
            stop("color_list must have at least k colors (one per cluster).");
        };
        cluster_colors <- stats::setNames(
            c(color_list[seq_len(k)], "grey70"),
            c(as.character(seq_len(k)), "Single")
        );
    };

    # Count per cluster and single
    cluster_counts <- table(data_plot$cluster[!duplicated(data_plot$Record_ID)]);

    cluster_labels <- stats::setNames(
        c(paste0("Cluster ", seq_len(k), " (", cluster_counts[as.character(seq_len(k))], ")"),
          paste0("Single point (", cluster_counts["Single"], ")")),
        c(as.character(seq_len(k)), "Single")
    );

    # ── 10. Helper: build trajectory plot with optional log x-axis ───────────
    .build_trajectory_plot <- function(log_x = FALSE) {
        data_single    <- data_plot %>% dplyr::filter(cluster == "Single");
        data_clustered <- data_plot %>% dplyr::filter(cluster != "Single");
        median_age <- median(data_plot$Age_VABS[data_plot$Instance == 1], na.rm = TRUE);

        p <- ggplot2::ggplot() +
            ggplot2::geom_point(
                data  = data_single,
                ggplot2::aes(x = Age_VABS, y = VABS_ABC, color = cluster),
                size  = 1.5,
                alpha = 0.6
            ) +
            ggplot2::geom_line(
                data      = data_clustered,
                ggplot2::aes(x = Age_VABS, y = VABS_ABC, group = Record_ID, color = cluster),
                alpha     = 0.4,
                linewidth = 0.45
            ) +
            ggplot2::geom_point(
                data  = data_clustered,
                ggplot2::aes(x = Age_VABS, y = VABS_ABC, color = cluster),
                size  = 1.2,
                alpha = 0.7
            ) +
            ggplot2::geom_vline(
                xintercept = median_age,
                color      = "red",
                linetype   = "dashed",
                linewidth  = 0.6
            ) +
            ggplot2::annotate(
                geom  = "text",
                x     = median_age,
                y     = max(data_plot$VABS_ABC, na.rm = TRUE),
                label = paste0("Median: ", round(median_age, 1)),
                color = "red",
                hjust = -0.1,
                vjust = 1.2,
                size  = 3.5
            ) +
            ggplot2::scale_color_manual(
                values = cluster_colors,
                labels = cluster_labels,
                breaks = c(as.character(seq_len(k)), "Single")
            ) +
            ggplot2::labs(
                x     = ifelse(log_x,
                               "Age at evaluation (years, log scale)",
                               "Age at evaluation (years)"),
                y     = "VABS ABC score",
                color = NULL
            ) +
            ggplot2::theme_classic(base_size = 13) +
            ggplot2::theme(
                legend.position      = c(0.98, 0.98),
                legend.justification = c("right", "top"),
                legend.background    = ggplot2::element_rect(fill = "white", color = "grey80", linewidth = 0.3),
                legend.margin        = ggplot2::margin(4, 6, 4, 6)
            );

        if (log_x) {
            p <- p + ggplot2::scale_x_log10();
        };

        return(p);
    };

    # ── 11. Build and save trajectory plots ───────────────────────────────────
    plot_trajectory        <- .build_trajectory_plot(log_x = FALSE);
    plot_trajectory_log    <- .build_trajectory_plot(log_x = TRUE);

    .save_plot(plot_trajectory,     "trajectory_by_cluster.png");
    .save_plot(plot_trajectory_log, "trajectory_by_cluster_logx.png");

    # ── 12. Cluster feature scatter plot ──────────────────────────────────────
    plot_cluster <- ggplot2::ggplot(
        cluster_data,
        ggplot2::aes(x = mean_vabs, y = slope_vabs, color = cluster)
    ) +
        ggplot2::geom_point(size = 2.5, alpha = 0.8) +
        ggplot2::scale_color_manual(
            values = cluster_colors[as.character(seq_len(k))],
            labels = cluster_labels[as.character(seq_len(k))]
        ) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        ggplot2::labs(
            title = "Clustering Features: Mean VABS vs. Slope",
            x     = "Mean VABS ABC (across all timepoints)",
            y     = "Slope of VABS ABC per year",
            color = "Cluster"
        ) +
        ggplot2::theme_bw(base_size = 13);

    .save_plot(plot_cluster, "cluster_features_scatter.png");

    # ── 13. Cluster summary statistics ────────────────────────────────────────
    cluster_summary <- cluster_data %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(
            n                  = dplyr::n(),
            mean_vabs_score    = round(mean(mean_vabs),     1),
            mean_baseline_vabs = round(mean(baseline_vabs), 1),
            mean_slope_vabs    = round(mean(slope_vabs),    3),
            mean_age_start     = round(mean(age_start),     1),
            .groups            = "drop"
        );

    # ── 14. Return all results (no plot objects) ───────────────────────────────
    return(list(
        model               = lme_model,
        data_with_residuals = df,
        patient_features    = cluster_data,
        data_plot           = data_plot,
        cluster_summary     = cluster_summary,
        n_clustered         = n_clustered,
        n_single            = n_single
    ));
}

VABS_trend.GBTM <- function(
        data,
        id_col        = "Record_ID",
        instance_col  = "Instance",
        age_col       = "Age_VABS",
        vabs_col      = "VABS_ABC",
        ng_max        = 4,
        log_age       = TRUE,
        color_list    = NULL,
        filepath      = NULL,
        plot_width    = 8,
        plot_height   = 6,
        plot_dpi      = 300,
        seed          = 42
    ) {

    # ── 0. Helper: save a ggplot ──────────────────────────────────────────────
    .save_plot <- function(plot_obj, filename) {
        if (!is.null(filepath)) {
            if (!dir.exists(filepath)) {
                dir.create(filepath, recursive = TRUE);
                message("Created directory: ", filepath);
            };
            full_path <- file.path(filepath, filename);
            ggplot2::ggsave(
                filename = full_path,
                plot     = plot_obj,
                width    = plot_width,
                height   = plot_height,
                dpi      = plot_dpi
            );
            message("Saved: ", full_path);
        };
    };

    # ── 1. Rename and filter ──────────────────────────────────────────────────
    df_all <- data %>%
        dplyr::rename(
            Record_ID = !!rlang::sym(id_col),
            Instance  = !!rlang::sym(instance_col),
            Age_VABS  = !!rlang::sym(age_col),
            VABS_ABC  = !!rlang::sym(vabs_col)
        ) %>%
        dplyr::filter(Instance > 0, !is.na(VABS_ABC), !is.na(Age_VABS)) %>%
        dplyr::arrange(Record_ID, Instance);

    # Exclude patients whose ALL VABS scores > 100
    ids_all_gt100 <- df_all %>%
        dplyr::group_by(Record_ID) %>%
        dplyr::filter(all(VABS_ABC > 100)) %>%
        dplyr::pull(Record_ID) %>%
        unique();

    message("Patients excluded (all VABS > 100): ", length(ids_all_gt100));

    df_all <- df_all %>% dplyr::filter(!Record_ID %in% ids_all_gt100);

    eligible_ids <- df_all %>%
        dplyr::count(Record_ID) %>%
        dplyr::filter(n >= 2) %>%
        dplyr::pull(Record_ID);

    df_model <- df_all %>% dplyr::filter(Record_ID %in% eligible_ids);

    message("Patients eligible for GBTM (>= 2 instances): ", length(eligible_ids));
    message("Patients with single instance (excluded):     ",
            dplyr::n_distinct(df_all$Record_ID) - length(eligible_ids));

    # ── 2. Integer subject ID and optional log-transform ─────────────────────
    df_model <- df_model %>%
        dplyr::mutate(Sample_ID = as.integer(factor(Record_ID)));

    df_model$Age_fit <- if (log_age) log(df_model$Age_VABS) else df_model$Age_VABS;

    # ── 3. Fit ng = 1 baseline model ─────────────────────────────────────────
    set.seed(seed);
    model1 <- lcmm::lcmm(
        VABS_ABC ~ Age_fit,
        random  = ~1,
        subject = "Sample_ID",
        ng      = 1,
        data    = df_model
    );

    # ── 4. Fit ng = 2 … ng_max models ────────────────────────────────────────
    models <- vector("list", ng_max);
    models[[1]] <- model1;

    for (g in seq(2, ng_max)) {
        set.seed(seed);
        models[[g]] <- lcmm::lcmm(
            VABS_ABC ~ Age_fit,
            random  = ~1,
            subject = "Sample_ID",
            ng      = g,
            mixture = ~Age_fit,
            nwg     = TRUE,
            data    = df_model,
            B       = model1
        );
    };

    # ── 5. Fit statistics ─────────────────────────────────────────────────────
    fit_stats <- purrr::map_dfr(seq_len(ng_max), function(g) {
        m <- models[[g]];
        entropy <- if (g == 1) {
            NA_real_;
        } else {
            probs  <- m$pprob[, grep("^prob", colnames(m$pprob)), drop = FALSE];
            n_subj <- nrow(probs);
            max_p  <- apply(probs, 1, max);
            H_max  <- log(g);
            1 - (-sum(max_p * log(pmax(max_p, 1e-10))) / (n_subj * H_max));
        };
        dplyr::tibble(
            G         = g,
            loglik    = m$loglik,
            AIC       = m$AIC,
            BIC       = m$BIC,
            entropy   = entropy,
            pct_class = if (g == 1) "100" else
                paste(round(100 * prop.table(table(m$pprob$class))), collapse = "/")
        );
    });

    message("\n── Model comparison ──────────────────────────────────────────");
    print(as.data.frame(fit_stats));

    best_g <- fit_stats %>%
        dplyr::filter(!is.na(entropy)) %>%
        dplyr::slice_max(entropy, n = 1, with_ties = FALSE) %>%
        dplyr::pull(G);

    message("\nBest model selected: ng = ", best_g, " (highest entropy)");

    # ── 6. Helper: assign classes and build data_plot for a given ng ──────────
    .make_data_plot <- function(g) {
        if (g == 1) {
            data_plot <- df_all %>%
                dplyr::left_join(
                    df_model %>% dplyr::distinct(Record_ID) %>%
                        dplyr::mutate(traj_class = "1"),
                    by = "Record_ID"
                ) %>%
                dplyr::mutate(
                    traj_class = dplyr::if_else(is.na(traj_class), "Single", traj_class),
                    traj_class = factor(traj_class, levels = c("Single", "1"))
                );
            return(data_plot);
        };

        m <- models[[g]];

        pprob_df <- m$pprob %>%
            dplyr::select(Sample_ID, class) %>%
            dplyr::rename(traj_class = class);

        id_map <- df_model %>% dplyr::distinct(Record_ID, Sample_ID);

        class_df <- id_map %>%
            dplyr::left_join(pprob_df, by = "Sample_ID") %>%
            dplyr::mutate(traj_class = as.character(traj_class));

        class_means <- df_model %>%
            dplyr::left_join(class_df, by = "Record_ID") %>%
            dplyr::group_by(traj_class) %>%
            dplyr::summarise(mean_vabs = mean(VABS_ABC), .groups = "drop") %>%
            dplyr::arrange(mean_vabs) %>%
            dplyr::mutate(new_label = as.character(dplyr::row_number()));

        recode_vec <- setNames(class_means$new_label, class_means$traj_class);
        class_df   <- class_df %>%
            dplyr::mutate(traj_class = recode_vec[traj_class]);

        data_plot <- df_all %>%
            dplyr::left_join(
                class_df %>% dplyr::select(Record_ID, traj_class),
                by = "Record_ID"
            ) %>%
            dplyr::mutate(
                traj_class = dplyr::if_else(is.na(traj_class), "Single", traj_class),
                traj_class = factor(traj_class,
                                    levels = c("Single", as.character(seq_len(g))))
            );
        return(data_plot);
    };

    # ── 7. Helper: build and save trajectory plots for a given ng ────────────
    .plot_for_g <- function(g) {
        data_plot <- .make_data_plot(g);

        if (is.null(color_list)) {
            traj_colors <- stats::setNames(
                c(grDevices::rainbow(g), "grey70"),
                c(as.character(seq_len(g)), "Single")
            );
        } else {
            if (length(color_list) < g) {
                warning("color_list has fewer than ", g,
                        " colors for ng=", g, "; falling back to rainbow.");
                traj_colors <- stats::setNames(
                    c(grDevices::rainbow(g), "grey70"),
                    c(as.character(seq_len(g)), "Single")
                );
            } else {
                traj_colors <- stats::setNames(
                    c(color_list[seq_len(g)], "grey70"),
                    c(as.character(seq_len(g)), "Single")
                );
            };
        };

        class_counts <- table(
            data_plot$traj_class[!duplicated(data_plot$Record_ID)]
        );

        traj_labels <- stats::setNames(
            c(paste0("Group ", seq_len(g),
                     " (", class_counts[as.character(seq_len(g))], ")"),
              paste0("Single point (", class_counts["Single"], ")")),
            c(as.character(seq_len(g)), "Single")
        );

        .build_one_plot <- function(log_x = FALSE) {
            data_single    <- data_plot %>% dplyr::filter(traj_class == "Single");
            data_clustered <- data_plot %>% dplyr::filter(traj_class != "Single");

            data_arrows <- data_clustered %>%
                dplyr::arrange(Record_ID, Age_VABS) %>%
                dplyr::group_by(Record_ID) %>%
                dplyr::mutate(
                    x_end = dplyr::lead(Age_VABS),
                    y_end = dplyr::lead(VABS_ABC)
                ) %>%
                dplyr::filter(!is.na(x_end)) %>%
                dplyr::ungroup();

            legend_pos     <- if (log_x) c(0.02, 0.98) else c(0.98, 0.98);
            legend_justify <- if (log_x) c("left", "top") else c("right", "top");

            p <- ggplot2::ggplot() +
                ggplot2::geom_point(
                    data  = data_single,
                    ggplot2::aes(x = Age_VABS, y = VABS_ABC, color = traj_class),
                    size  = 1.5
                ) +
                ggplot2::geom_segment(
                    data = data_arrows,
                    ggplot2::aes(x = Age_VABS, y = VABS_ABC,
                                 xend = x_end, yend = y_end,
                                 color = traj_class),
                    arrow     = ggplot2::arrow(length = ggplot2::unit(0.15, "cm"),
                                               type = "closed"),
                    linewidth = 0.45
                ) +
                ggplot2::geom_point(
                    data  = data_clustered,
                    ggplot2::aes(x = Age_VABS, y = VABS_ABC, color = traj_class),
                    size  = 1.2
                ) +
                ggplot2::geom_vline(
                    xintercept = 10,
                    color      = "red",
                    linetype   = "dashed",
                    linewidth  = 0.6
                ) +
                ggplot2::scale_color_manual(
                    values = traj_colors,
                    labels = traj_labels,
                    breaks = c(as.character(seq_len(g)), "Single")
                ) +
                ggplot2::labs(
                    title = NULL,
                    x     = ifelse(log_x,
                                   "Age at evaluation (years, log scale)",
                                   "Age at evaluation (years)"),
                    y     = "VABS ABC score",
                    color = NULL
                ) +
                ggplot2::theme_classic(base_size = 13) +
                ggplot2::theme(
                    legend.position      = legend_pos,
                    legend.justification = legend_justify,
                    legend.background    = ggplot2::element_blank(),
                    legend.key           = ggplot2::element_blank(),
                    legend.margin        = ggplot2::margin(4, 6, 4, 6)
                );

            if (log_x) {
                p <- p + ggplot2::scale_x_log10();
            };
            return(p);
        };

        p_linear <- .build_one_plot(log_x = FALSE);
        p_logx   <- .build_one_plot(log_x = TRUE);

        tag <- if (g == best_g) paste0("ng", g, "_BEST") else paste0("ng", g);
        .save_plot(p_linear, paste0("trajectory_GBTM_", tag, ".png"));
        .save_plot(p_logx,   paste0("trajectory_GBTM_", tag, "_logx.png"));
    };

    # ── 8. Generate plots for ng = 2 … ng_max ────────────────────────────────
    for (g in seq(2, ng_max)) {
        message("Plotting ng = ", g, " ...");
        .plot_for_g(g);
    };

    # ── 9. Append Record_ID to Sample_ID-based model outputs ─────────────────
    id_map <- df_model %>% dplyr::distinct(Record_ID, Sample_ID);

    for (g in seq_len(ng_max)) {
        for (slot in c("pred", "pprob", "predRE", "classpredRE")) {
            if (!is.null(models[[g]][[slot]])) {
                models[[g]][[slot]] <- models[[g]][[slot]] %>%
                    dplyr::left_join(id_map, by = "Sample_ID") %>%
                    dplyr::relocate(Record_ID, .after = Sample_ID);
            };
        };
    };

    # ── 10. Return ────────────────────────────────────────────────────────────
    return(models);
}

molecular_data.heatmap <- function(molecular_data, clinical_data, type = "h",
                                   height = 8, width = 16, filename = NULL) {

    # Validate type parameter: "h" = horizontal, "v" = vertical
    if (!type %in% c("h", "v")) {
        stop("type must be 'h' (horizontal) or 'v' (vertical)");
    };

    is_horizontal <- type == "h";

    # Transpose the matrix: rows = molecular features, columns = Variant_name
    mat <- t(molecular_data);
    variant_names <- colnames(mat);

    # Collect VABS_ABC values for each Variant_name (multiple patients may share a variant)
    vabs_list <- lapply(variant_names, function(vn) {
        vals <- clinical_data$VABS_ABC[clinical_data$Variant_name == vn];
        vals <- vals[!is.na(vals)];
        if (length(vals) == 0) return(NULL);
        vals;
    });
    names(vabs_list) <- variant_names;

    # Compute global VABS_ABC range for a unified axis across all boxplots
    all_vabs <- unlist(vabs_list);
    vabs_range <- range(all_vabs, na.rm = TRUE);
    vabs_ylim <- c(floor(vabs_range[1]) - 5, ceiling(vabs_range[2]) + 5);

    # Compute pretty tick positions for gridlines
    vabs_ticks <- pretty(vabs_ylim, n = 5);
    vabs_ticks <- vabs_ticks[vabs_ticks >= vabs_ylim[1] & vabs_ticks <= vabs_ylim[2]];

    # Alternating background colors for each variant strip
    bg_colors <- c("white", "#F0F0F0");

    # Boxplot style constants
    line_col    <- "#1F77B4";   # Line and border color
    fill_col    <- "#AEC7E8";   # Box fill color
    box_hw      <- 0.35;        # Half-width (h) or half-height (v) of box
    whisker_lwd <- 3.0;         # Whisker line width (bold)
    border_lwd  <- whisker_lwd; # Box border matches whisker line width
    iqr_lwd     <- 1.5;         # Q1/Q3 edge line width
    median_lwd  <- 3.0;         # Median line width (bold black)

    # Cluster colors for 4 groups
    cluster_colors <- c("#377EB8", "#4DAF4A", "#FF7F00", "#E41A1C");

    # Helper: reorder dendrogram to achieve target cluster order
    rotate_dend <- function(data_mat) {
        hc   <- hclust(dist(data_mat), method = "ward.D2");
        dend <- as.dendrogram(hc);

        # Step 1: swap children within the left subtree of root
        dend[[1]] <- dendextend::rotate(dend[[1]], order = rev(labels(dend[[1]])));

        # Step 2: swap the two children of the root node
        dend <- dendextend::rotate(dend, order = rev(labels(dend)));

        # Convert rotated dendrogram back to hclust
        hc_rotated <- as.hclust(dend);

        list(dend = dend, hclust = hc_rotated);
    };

    # Unified annotation drawing function used by AnnotationFunction
    # draw_one_boxplot is defined inside anno_draw to ensure it is accessible
    # within the AnnotationFunction execution environment
    anno_draw <- function(index, ...) {

        # In vertical left_annotation mode, index runs bottom-to-top;
        # reverse it so boxplots align with rows top-to-bottom
        if (!is_horizontal) {
            index <- rev(index);
        };

        # Draw one boxplot for a single variant
        # h mode: index i maps to x position, VABS on y-axis
        # v mode: index i maps to y position, VABS on x-axis
        draw_one_boxplot <- function(i, vals) {

            if (length(vals) < 3) {
                # Fewer than 3 patients: draw individual solid points only
                if (is_horizontal) {
                    grid.points(x = rep(i, length(vals)), y = vals,
                                pch = 16, size = unit(2, "mm"),
                                gp = gpar(col = "black"), default.units = "native");
                } else {
                    grid.points(x = vals, y = rep(i, length(vals)),
                                pch = 16, size = unit(2, "mm"),
                                gp = gpar(col = "black"), default.units = "native");
                };
                return(invisible(NULL));
            };

            # Compute boxplot statistics
            q          <- quantile(vals, probs = c(0.25, 0.5, 0.75));
            iqr        <- q[3] - q[1];
            whisker_lo <- max(min(vals), q[1] - 1.5 * iqr);
            whisker_hi <- min(max(vals), q[3] + 1.5 * iqr);
            outliers   <- vals[vals < whisker_lo | vals > whisker_hi];

            if (is_horizontal) {
                # IQR box
                grid.rect(x = unit(i, "native"), y = unit(q[1], "native"),
                          width = unit(box_hw * 2, "native"), height = unit(q[3] - q[1], "native"),
                          just = c("center", "bottom"),
                          gp = gpar(fill = fill_col, col = line_col, lwd = border_lwd));
                # Q1 and Q3 horizontal edges
                grid.segments(x0 = unit(i - box_hw, "native"), x1 = unit(i + box_hw, "native"),
                              y0 = unit(q[1], "native"),        y1 = unit(q[1], "native"),
                              gp = gpar(col = line_col, lwd = iqr_lwd));
                grid.segments(x0 = unit(i - box_hw, "native"), x1 = unit(i + box_hw, "native"),
                              y0 = unit(q[3], "native"),        y1 = unit(q[3], "native"),
                              gp = gpar(col = line_col, lwd = iqr_lwd));
                # Median line: bold black horizontal
                grid.segments(x0 = unit(i - box_hw, "native"), x1 = unit(i + box_hw, "native"),
                              y0 = unit(q[2], "native"),        y1 = unit(q[2], "native"),
                              gp = gpar(col = "black", lwd = median_lwd));
                # Lower and upper whiskers: vertical lines
                grid.segments(x0 = unit(i, "native"), x1 = unit(i, "native"),
                              y0 = unit(whisker_lo, "native"), y1 = unit(q[1], "native"),
                              gp = gpar(col = line_col, lwd = whisker_lwd));
                grid.segments(x0 = unit(i, "native"), x1 = unit(i, "native"),
                              y0 = unit(q[3], "native"), y1 = unit(whisker_hi, "native"),
                              gp = gpar(col = line_col, lwd = whisker_lwd));
                # Outliers
                if (length(outliers) > 0) {
                    grid.points(x = rep(i, length(outliers)), y = outliers,
                                pch = 16, size = unit(1.5, "mm"),
                                gp = gpar(col = "black"), default.units = "native");
                };

            } else {
                # IQR box (drawn horizontally)
                grid.rect(x = unit(q[1], "native"), y = unit(i, "native"),
                          width = unit(q[3] - q[1], "native"), height = unit(box_hw * 2, "native"),
                          just = c("left", "center"),
                          gp = gpar(fill = fill_col, col = line_col, lwd = border_lwd));
                # Q1 and Q3 vertical edges
                grid.segments(x0 = unit(q[1], "native"), x1 = unit(q[1], "native"),
                              y0 = unit(i - box_hw, "native"), y1 = unit(i + box_hw, "native"),
                              gp = gpar(col = line_col, lwd = iqr_lwd));
                grid.segments(x0 = unit(q[3], "native"), x1 = unit(q[3], "native"),
                              y0 = unit(i - box_hw, "native"), y1 = unit(i + box_hw, "native"),
                              gp = gpar(col = line_col, lwd = iqr_lwd));
                # Median line: bold black vertical
                grid.segments(x0 = unit(q[2], "native"), x1 = unit(q[2], "native"),
                              y0 = unit(i - box_hw, "native"), y1 = unit(i + box_hw, "native"),
                              gp = gpar(col = "black", lwd = median_lwd));
                # Left and right whiskers: horizontal lines
                grid.segments(x0 = unit(whisker_lo, "native"), x1 = unit(q[1], "native"),
                              y0 = unit(i, "native"),           y1 = unit(i, "native"),
                              gp = gpar(col = line_col, lwd = whisker_lwd));
                grid.segments(x0 = unit(q[3], "native"),       x1 = unit(whisker_hi, "native"),
                              y0 = unit(i, "native"),           y1 = unit(i, "native"),
                              gp = gpar(col = line_col, lwd = whisker_lwd));
                # Outliers
                if (length(outliers) > 0) {
                    grid.points(x = outliers, y = rep(i, length(outliers)),
                                pch = 16, size = unit(1.5, "mm"),
                                gp = gpar(col = "black"), default.units = "native");
                };
            };
        };

        n <- length(index);

        # Set up viewport: VABS axis vs strip index axis
        if (is_horizontal) {
            pushViewport(viewport(xscale = c(0.5, n + 0.5), yscale = vabs_ylim));
        } else {
            pushViewport(viewport(xscale = vabs_ylim, yscale = c(0.5, n + 0.5)));
        };

        # Draw alternating strip backgrounds
        for (i in seq_len(n)) {
            bg <- gpar(fill = bg_colors[(i %% 2) + 1], col = NA);
            if (is_horizontal) {
                grid.rect(x = unit(i, "native"), y = unit(vabs_ylim[1], "native"),
                          width = unit(1, "native"), height = unit(diff(vabs_ylim), "native"),
                          just = c("center", "bottom"), gp = bg);
            } else {
                grid.rect(x = unit(vabs_ylim[1], "native"), y = unit(i, "native"),
                          width = unit(diff(vabs_ylim), "native"), height = unit(1, "native"),
                          just = c("left", "center"), gp = bg);
            };
        };

        # Draw dashed gridlines along the VABS axis
        for (tick in vabs_ticks) {
            if (is_horizontal) {
                grid.segments(x0 = unit(0.5, "native"),  x1 = unit(n + 0.5, "native"),
                              y0 = unit(tick, "native"), y1 = unit(tick, "native"),
                              gp = gpar(col = "gray70", lty = "dashed", lwd = 0.6));
            } else {
                grid.segments(x0 = unit(tick, "native"), x1 = unit(tick, "native"),
                              y0 = unit(0.5, "native"),  y1 = unit(n + 0.5, "native"),
                              gp = gpar(col = "gray70", lty = "dashed", lwd = 0.6));
            };
        };

        # Draw boxplot for each variant in index
        for (i in seq_along(index)) {
            vn   <- variant_names[index[i]];
            vals <- vabs_list[[vn]];
            if (!is.null(vals) && length(vals) > 0) {
                draw_one_boxplot(i, vals);
            };
            # If vabs_list[[vn]] is NULL, leave the strip blank
        };

        # Draw axis and border
        if (is_horizontal) {
            grid.yaxis(at = vabs_ticks, gp = gpar(fontsize = 7));
        } else {
            grid.xaxis(at = vabs_ticks, gp = gpar(fontsize = 7), 
                       edits = gEdit("labels", y = unit(-1.5, "mm"), just = "top"));
        };
        grid.rect(gp = gpar(col = "black", fill = NA, lwd = 1));
        upViewport();
    };

    # Build annotation object depending on type
    if (is_horizontal) {
        anno <- HeatmapAnnotation(
            VABS_ABC = AnnotationFunction(
                fun    = anno_draw,
                height = unit(4, "cm"),
                n      = ncol(mat)
            ),
            annotation_label       = "Final VABS ABC",
            annotation_name_side   = "left",
            annotation_name_rot    = 90,
            annotation_name_gp     = gpar(fontsize = 9),
            annotation_name_offset = unit(8, "mm")
        );
    } else {
        anno <- rowAnnotation(
            VABS_ABC = AnnotationFunction(
                fun   = anno_draw,
                width = unit(4, "cm"),
                n     = nrow(t(mat)),
                which = "row"
            ),
            space = anno_empty(border = FALSE, width = unit(0.1, "cm")), 
            annotation_label       = "Final VABS ABC",
            annotation_name_side   = "bottom",
            annotation_name_rot    = 0,
            annotation_name_gp     = gpar(fontsize = 9),
            annotation_name_offset = unit(5, "mm")
        );
    };

    # Color mapping: reversed RdYlBu palette (blue = low, yellow = mid, red = high)
    col_fun <- colorRamp2(
        seq(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE), length.out = 11),
        rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu"))
    );

    if (is_horizontal) {
        # Compute rotated dendrogram and hclust for columns (variants)
        col_result <- rotate_dend(t(mat));
        col_dend   <- dendextend::color_branches(col_result$dend, k = 4, col = cluster_colors);

        ht <- Heatmap(
            mat,
            name                 = "Scaled\nValue",
            col                  = col_fun,
            top_annotation       = anno,
            cluster_rows         = FALSE,
            cluster_columns      = col_dend,
            show_column_names    = TRUE,
            show_row_names       = TRUE,
            column_names_gp      = gpar(fontsize = 8),
            row_names_gp         = gpar(fontsize = 9),
            column_names_rot     = 45,
            rect_gp              = gpar(col = "gray60", lwd = 0.5),
            heatmap_legend_param = list(
                title     = "Scaled\nValue",
                title_gp  = gpar(fontsize = 9),
                labels_gp = gpar(fontsize = 8)
            )
        );

    } else {
        # Compute rotated dendrogram and hclust for rows (variants)
        mat_v      <- t(mat);
        row_result <- rotate_dend(mat_v);
        row_dend   <- dendextend::color_branches(row_result$dend, k = 4, col = cluster_colors);

        ht <- Heatmap(
            mat_v,
            name                 = "Scaled\nValue",
            col                  = col_fun,
            left_annotation      = anno,
            cluster_rows         = row_dend,
            cluster_columns      = FALSE,
            show_column_names    = TRUE,
            show_row_names       = TRUE,
            column_names_gp      = gpar(fontsize = 8),
            row_names_gp         = gpar(fontsize = 9),
            column_names_rot     = 45,
            rect_gp              = gpar(col = "gray60", lwd = 0.5),
            heatmap_legend_param = list(
                title     = "Scaled\nValue",
                title_gp  = gpar(fontsize = 9),
                labels_gp = gpar(fontsize = 8)
            )
        );
    };

    # Draw and save or return
    if (!is.null(filename)) {
        ext <- tolower(tools::file_ext(filename));
        if (ext == "pdf") {
            pdf(filename, width = width, height = height);
        } else if (ext == "svg") {
            svg(filename, width = width, height = height);
        } else if (ext == "png") {
            png(filename, width = width, height = height, units = "in", res = 300);
        } else if (ext %in% c("tif", "tiff")) {
            tiff(filename, width = width, height = height, units = "in", res = 300);
        } else {
            stop(paste("Unsupported file extension:", ext));
        };
        draw(ht, heatmap_legend_side = "right");
        dev.off();
        message("Heatmap saved to: ", filename);
    } else {
        draw(ht, heatmap_legend_side = "right");
    };

    # Always return ht and rotated hclust invisibly
    if (is_horizontal) {
        invisible(list(ht = ht, hclust = col_result$hclust));
    } else {
        invisible(list(ht = ht, hclust = row_result$hclust));
    };
}

regression_model <- function(
    data,                        # Data frame containing all variables
    predictor,                   # Character vector of predictor variable names
    response,                    # Character string of response variable name
    model        = c("en", "rf", "xgb"),  # Model type: "en" (Elastic Net), "rf" (Random Forest), "xgb" (XGBoost)
    id_col       = "Record_ID",  # Name of the sample ID column
    variant_col  = "Variant_name", # Name of the variant column for grouping
    filepath     = NULL,         # Directory path to save all outputs (default: working directory)
    model_name   = NULL,         # Prefix for all saved output files (default: auto-generated from model type)
    tune_folds   = 5,            # Number of CV folds for global parameter tuning
    tune_repeats = 3,            # Number of repeats for repeated CV tuning
    seed         = 123,          # Random seed for reproducibility
    n_cores      = 1,            # Number of parallel cores for tuning (>1 requires doParallel)
    save_plots   = TRUE,         # Whether to save diagnostic plots as PNG files
    verbose      = TRUE          # Whether to print progress messages
    ) {
    # ── Validate arguments ──────────────────────────────────────────────────
    model <- match.arg(model);
    if (is.null(model_name)) model_name <- paste0("lovo_", model);
    if (is.null(filepath))   filepath   <- getwd();
    if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE);

    msg <- function(...) if (verbose) message(...);

    # ── Load required packages ──────────────────────────────────────────────
    required_pkgs <- c("caret", "ggplot2", "dplyr");
    model_pkgs <- switch(model,
        en  = c("glmnet"),
        rf  = c("randomForest"),
        xgb = c("xgboost", "Matrix")
    );
    for (pkg in c(required_pkgs, model_pkgs)) {
        if (!requireNamespace(pkg, quietly = TRUE))
            stop("Package '", pkg, "' is required but not installed.");
    }
    library(caret, quietly = TRUE);

    # ── Parallel backend ────────────────────────────────────────────────────
    if (n_cores > 1) {
        if (!requireNamespace("doParallel", quietly = TRUE))
            stop("Package 'doParallel' is required for parallel computing.");
        cl <- parallel::makeCluster(n_cores);
        doParallel::registerDoParallel(cl);
        on.exit(parallel::stopCluster(cl), add = TRUE);
        msg("  Using ", n_cores, " parallel cores.");
    }

    # ── Prepare data ────────────────────────────────────────────────────────
    msg("=", strrep("=", 59));
    msg("  Model : ", toupper(model));
    msg("=", strrep("=", 59));
    msg("\n[1] Preparing data ...");

    # Check variant column exists
    if (!variant_col %in% colnames(data))
        stop("Column '", variant_col, "' not found in data.");

    # Keep only required columns; drop rows with any NA in predictor/response
    keep_cols <- c(id_col, variant_col, predictor, response);
    temp      <- data[, intersect(keep_cols, colnames(data)), drop = FALSE];
    complete  <- complete.cases(temp[, c(predictor, response)]);
    temp      <- temp[complete, ];
    n         <- nrow(temp);
    msg("    Samples after NA removal : ", n);
    msg("    Predictors               : ", length(predictor));

    # Get unique variants and their sample indices
    unique_variants <- unique(temp[[variant_col]]);
    n_variants      <- length(unique_variants);
    msg("    Unique variants          : ", n_variants);

    # Print variant size distribution
    variant_sizes <- table(temp[[variant_col]]);
    msg("    Samples per variant      : min=", min(variant_sizes),
        ", median=", median(variant_sizes),
        ", max=", max(variant_sizes));

    sample_ids <- temp[[id_col]];
    X_raw      <- temp[, predictor, drop = FALSE];
    y          <- temp[[response]];

    # ── Encode categorical variables ─────────────────────────────────────────
    cat_cols <- predictor[sapply(X_raw, function(v) is.factor(v) || is.character(v))];
    num_cols <- setdiff(predictor, cat_cols);

    if (length(cat_cols) > 0) {
        msg("    Categorical predictors   : ", paste(cat_cols, collapse = ", "));
        dummy_formula <- as.formula(paste("~", paste(cat_cols, collapse = " + ")));
        dummy_mat     <- model.matrix(dummy_formula, data = X_raw)[, -1, drop = FALSE];
        X_proc        <- cbind(X_raw[, num_cols, drop = FALSE], dummy_mat);
    } else {
        X_proc <- X_raw;
    }
    X_proc <- as.data.frame(lapply(X_proc, as.numeric));
    msg("    Final feature count      : ", ncol(X_proc));

    # ── Step 1: Global parameter tuning ─────────────────────────────────────
    msg("\n[2] Global parameter tuning (", tune_folds, "-fold CV x", tune_repeats, " repeats) ...");
    set.seed(seed);
    train_ctrl <- trainControl(
        method        = "repeatedcv",
        number        = tune_folds,
        repeats       = tune_repeats,
        allowParallel = (n_cores > 1)
    );

    if (model == "en") {
        tune_grid <- expand.grid(
            alpha  = seq(0, 1, by = 0.2),
            lambda = 10^seq(-4, 1, length.out = 20)
        );
        set.seed(seed);
        tuned <- caret::train(
            x          = X_proc,
            y          = y,
            method     = "glmnet",
            trControl  = train_ctrl,
            tuneGrid   = tune_grid,
            preProcess = c("center", "scale")
        );

    } else if (model == "rf") {
        tune_grid <- data.frame(
            mtry = unique(round(seq(2, max(2, ncol(X_proc)), length.out = 8)))
        );
        set.seed(seed);
        tuned <- caret::train(
            x          = X_proc,
            y          = y,
            method     = "rf",
            trControl  = train_ctrl,
            tuneGrid   = tune_grid,
            ntree      = 500,
            importance = TRUE
        );

    } else if (model == "xgb") {
        tune_grid <- expand.grid(
            nrounds          = c(100, 200, 300),
            max_depth        = c(3, 5, 7),
            eta              = c(0.01, 0.05, 0.1),
            gamma            = 0,
            colsample_bytree = c(0.7, 1.0),
            min_child_weight = 1,
            subsample        = 0.8
        );
        set.seed(seed);
        tuned <- caret::train(
            x         = X_proc,
            y         = y,
            method    = "xgbTree",
            trControl = train_ctrl,
            tuneGrid  = tune_grid,
            verbosity = 0
        );
    }

    best_params <- tuned$bestTune;
    msg("    Best parameters:");
    for (nm in names(best_params))
        msg("      ", nm, " = ", best_params[[nm]]);

    # Save tuning plot
    if (save_plots) {
        tune_plot_file <- file.path(filepath, paste0(model_name, ".tuning.png"));
        png(tune_plot_file, width = 7, height = 5, units = "in", res = 300);
        print(plot(tuned, main = paste(toupper(model), "- Parameter Tuning")));
        dev.off();
        msg("    Tuning plot saved : ", tune_plot_file);
    }

    # ── Step 2: Leave-One-Variant-Out CV (LOVO-CV) ──────────────────────────
    msg("\n[3] Running Leave-One-Variant-Out CV (", n_variants, " variants, ",
        n, " samples) with fixed best parameters ...");
    msg("    Each iteration withholds ALL samples carrying one unique variant.");

    lovo_pred <- data.frame(
        SampleID         = sample_ids,
        Variant          = temp[[variant_col]],
        Actual           = y,
        Predicted        = NA_real_,
        stringsAsFactors = FALSE
    );

    # Helper: train on training samples, predict test samples
    .fit_predict <- function(train_X, train_y, test_X) {
        if (model == "en") {
            library(glmnet, quietly = TRUE);
            pre  <- caret::preProcess(train_X, method = c("center", "scale"));
            trX  <- predict(pre, train_X);
            teX  <- predict(pre, test_X);
            fit  <- glmnet::glmnet(
                x      = as.matrix(trX),
                y      = train_y,
                alpha  = best_params$alpha,
                lambda = best_params$lambda
            );
            pred <- as.numeric(predict(fit, newx = as.matrix(teX), s = best_params$lambda));

        } else if (model == "rf") {
            library(randomForest, quietly = TRUE);
            fit  <- randomForest::randomForest(
                x          = train_X,
                y          = train_y,
                ntree      = 500,
                mtry       = best_params$mtry,
                importance = FALSE
            );
            pred <- as.numeric(predict(fit, newdata = test_X));

        } else if (model == "xgb") {
            library(xgboost, quietly = TRUE);
            dtrain <- xgboost::xgb.DMatrix(
                data  = as.matrix(train_X),
                label = train_y
            );
            dtest <- xgboost::xgb.DMatrix(data = as.matrix(test_X));
            fit   <- xgboost::xgb.train(
                params = list(
                    objective        = "reg:squarederror",
                    max_depth        = best_params$max_depth,
                    eta              = best_params$eta,
                    gamma            = best_params$gamma,
                    colsample_bytree = best_params$colsample_bytree,
                    min_child_weight = best_params$min_child_weight,
                    subsample        = best_params$subsample,
                    verbosity        = 0
                ),
                data    = dtrain,
                nrounds = best_params$nrounds
            );
            pred <- as.numeric(predict(fit, dtest));
        }
        return(pred);
    }

    # ── LOVO-CV loop: iterate over unique variants ───────────────────────────
    for (v_idx in seq_len(n_variants)) {
        variant_name <- unique_variants[v_idx];

        # Boolean index: which rows carry this variant
        test_mask  <- temp[[variant_col]] == variant_name;
        train_mask <- !test_mask;

        n_test  <- sum(test_mask);
        n_train <- sum(train_mask);

        if (verbose && (v_idx %% 20 == 0 || v_idx == 1 || v_idx == n_variants))
            msg("    Variant ", v_idx, " / ", n_variants,
                " : '", variant_name, "'",
                "  (train=", n_train, ", test=", n_test, ")");

        # Safety check: skip if training set is too small
        if (n_train < 5) {
            warning("Variant '", variant_name, "' skipped: only ", n_train,
                    " training samples remaining. Predictions set to NA.");
            next;
        }

        train_X <- X_proc[train_mask, , drop = FALSE];
        train_y <- y[train_mask];
        test_X  <- X_proc[test_mask,  , drop = FALSE];

        set.seed(seed + v_idx);
        preds <- .fit_predict(train_X, train_y, test_X);

        # Assign predictions back to all rows carrying this variant
        lovo_pred$Predicted[test_mask] <- preds;
    }

    # ── Step 3: Compute performance metrics ──────────────────────────────────
    msg("\n[4] Computing performance metrics ...");

    # Only use samples that received a prediction (non-NA)
    valid_idx <- !is.na(lovo_pred$Predicted);
    n_valid   <- sum(valid_idx);
    msg("    Samples with predictions : ", n_valid, " / ", n);

    resid   <- lovo_pred$Actual[valid_idx] - lovo_pred$Predicted[valid_idx];
    rmse    <- sqrt(mean(resid^2));
    mae     <- mean(abs(resid));
    ss_res  <- sum(resid^2);
    ss_tot  <- sum((lovo_pred$Actual[valid_idx] - mean(lovo_pred$Actual[valid_idx]))^2);
    r2      <- 1 - ss_res / ss_tot;
    pearson <- cor(lovo_pred$Actual[valid_idx], lovo_pred$Predicted[valid_idx],
                   method = "pearson");

    metrics <- data.frame(
        Model    = model,
        CV_type  = "LOVO",
        N        = n_valid,
        Variants = n_variants,
        RMSE     = round(rmse,    4),
        MAE      = round(mae,     4),
        R2       = round(r2,      4),
        Pearson  = round(pearson, 4)
    );
    msg("    RMSE    = ", round(rmse,    4));
    msg("    MAE     = ", round(mae,     4));
    msg("    R2      = ", round(r2,      4));
    msg("    Pearson = ", round(pearson, 4));

    # ── Step 4: Feature importance on full dataset ───────────────────────────
    msg("\n[5] Computing feature importance on full dataset ...");
    feature_importance <- NULL;

    if (model == "en") {
        library(glmnet, quietly = TRUE);
        pre_full <- caret::preProcess(X_proc, method = c("center", "scale"));
        X_scaled <- predict(pre_full, X_proc);
        fit_full <- glmnet::glmnet(
            x      = as.matrix(X_scaled),
            y      = y,
            alpha  = best_params$alpha,
            lambda = best_params$lambda
        );
        coefs <- as.matrix(coef(fit_full, s = best_params$lambda));
        feature_importance <- data.frame(
            Feature     = rownames(coefs)[-1],
            Importance  = abs(coefs[-1, 1]),
            Coefficient = coefs[-1, 1]
        );
        feature_importance <- feature_importance[
            order(feature_importance$Importance, decreasing = TRUE), ];

    } else if (model == "rf") {
        library(randomForest, quietly = TRUE);
        fit_full <- randomForest::randomForest(
            x          = X_proc,
            y          = y,
            ntree      = 500,
            mtry       = best_params$mtry,
            importance = TRUE
        );
        imp_mat <- randomForest::importance(fit_full);
        feature_importance <- data.frame(
            Feature       = rownames(imp_mat),
            Importance    = imp_mat[, "%IncMSE"],
            IncNodePurity = imp_mat[, "IncNodePurity"]
        );
        feature_importance <- feature_importance[
            order(feature_importance$Importance, decreasing = TRUE), ];

    } else if (model == "xgb") {
        library(xgboost, quietly = TRUE);
        dtrain_full <- xgboost::xgb.DMatrix(
            data  = as.matrix(X_proc),
            label = y
        );
        fit_full <- xgboost::xgb.train(
            params = list(
                objective        = "reg:squarederror",
                max_depth        = best_params$max_depth,
                eta              = best_params$eta,
                gamma            = best_params$gamma,
                colsample_bytree = best_params$colsample_bytree,
                min_child_weight = best_params$min_child_weight,
                subsample        = best_params$subsample,
                verbosity        = 0
            ),
            data    = dtrain_full,
            nrounds = best_params$nrounds
        );
        imp_raw <- xgboost::xgb.importance(
            feature_names = colnames(X_proc),
            model         = fit_full
        );
        feature_importance <- as.data.frame(imp_raw);
        colnames(feature_importance)[1] <- "Feature";
        colnames(feature_importance)[2] <- "Importance";
    }

    # ── Step 5: Save outputs ─────────────────────────────────────────────────
    msg("\n[6] Saving outputs to : ", filepath);

    # 5a. LOVO predictions (includes Variant column)
    pred_file <- file.path(filepath, paste0(model_name, ".lovo_predictions.txt"));
    write.table(lovo_pred, file = pred_file, row.names = FALSE, quote = FALSE, sep = "\t");
    msg("    Predictions saved  : ", pred_file);

    # 5b. Performance metrics
    metric_file <- file.path(filepath, paste0(model_name, ".metrics.txt"));
    write.table(metrics, file = metric_file, row.names = FALSE, quote = FALSE, sep = "\t");
    msg("    Metrics saved      : ", metric_file);

    # 5c. Feature importance table
    if (!is.null(feature_importance)) {
        imp_file <- file.path(filepath, paste0(model_name, ".feature_importance.txt"));
        write.table(feature_importance, file = imp_file, row.names = FALSE, quote = FALSE, sep = "\t");
        msg("    Importance saved   : ", imp_file);
    }

    # 5d. Diagnostic plots
    if (save_plots) {
        valid_pred <- lovo_pred[valid_idx, ];

        # (i) Actual vs Predicted scatter, colored by variant group size
        ap_file <- file.path(filepath, paste0(model_name, ".actual_vs_predicted.png"));
        png(ap_file, width = 5, height = 5, units = "in", res = 300);
        plot(
            valid_pred$Actual, valid_pred$Predicted,
            xlab = "Actual", ylab = "Predicted",
            main = paste(toupper(model), "- LOVO-CV: Actual vs Predicted"),
            pch  = 19, col  = "#2166AC88"
        );
        abline(a = 0, b = 1, col = "red", lty = 2, lwd = 1.5);
        legend("topleft", bty = "n",
               legend = c(
                   paste0("RMSE = ", round(rmse,    3)),
                   paste0("R\u00b2 = ",   round(r2,      3)),
                   paste0("r = ",    round(pearson, 3)),
                   paste0("N variants = ", n_variants)
               ));
        dev.off();
        msg("    Actual vs predicted plot saved : ", ap_file);

        # (ii) Residual plot
        res_file <- file.path(filepath, paste0(model_name, ".residuals.png"));
        resid_valid <- valid_pred$Actual - valid_pred$Predicted;
        png(res_file, width = 5, height = 4, units = "in", res = 300);
        plot(
            valid_pred$Predicted, resid_valid,
            xlab = "Predicted", ylab = "Residual",
            main = paste(toupper(model), "- LOVO-CV Residuals"),
            pch  = 19, col  = "#D6604D88"
        );
        abline(h = 0, col = "gray40", lty = 2);
        dev.off();
        msg("    Residual plot saved            : ", res_file);

        # (iii) Per-variant prediction summary plot
        variant_summary <- do.call(rbind, lapply(unique_variants, function(vn) {
            rows <- valid_pred[valid_pred$Variant == vn, ];
            if (nrow(rows) == 0) return(NULL);
            data.frame(
                Variant    = vn,
                N          = nrow(rows),
                Mean_Actual    = mean(rows$Actual),
                Mean_Predicted = mean(rows$Predicted),
                stringsAsFactors = FALSE
            );
        }));
        variant_summary <- variant_summary[!is.null(variant_summary), ];

        vs_file <- file.path(filepath, paste0(model_name, ".variant_summary.png"));
        png(vs_file, width = 6, height = 5, units = "in", res = 300);
        plot(
            variant_summary$Mean_Actual, variant_summary$Mean_Predicted,
            xlab = "Mean Actual (per variant)",
            ylab = "Mean Predicted (per variant)",
            main = paste(toupper(model), "- LOVO-CV: Per-Variant Summary"),
            pch  = 21, bg = "#4393C388", col = "#2166AC",
            cex  = pmax(0.5, log1p(variant_summary$N) * 0.6)
        );
        abline(a = 0, b = 1, col = "red", lty = 2, lwd = 1.5);
        legend("topleft", bty = "n",
               legend = c("Point size ∝ log(N samples)", paste0("N variants = ", nrow(variant_summary))));
        dev.off();
        msg("    Variant summary plot saved     : ", vs_file);

        # (iv) Feature importance bar plot (top 20)
        if (!is.null(feature_importance) && nrow(feature_importance) > 0) {
            top_n  <- min(20, nrow(feature_importance));
            top_df <- feature_importance[seq_len(top_n), ];
            top_df$Feature <- factor(top_df$Feature, levels = rev(top_df$Feature));
            p_imp <- ggplot2::ggplot(top_df, ggplot2::aes(x = Feature, y = Importance)) +
                ggplot2::geom_bar(stat = "identity", fill = "#4393C3") +
                ggplot2::coord_flip() +
                ggplot2::labs(
                    title = paste(toupper(model), "- Feature Importance (Top", top_n, ")"),
                    x = NULL, y = "Importance"
                ) +
                ggplot2::theme_bw(base_size = 11);
            imp_plot_file <- file.path(filepath, paste0(model_name, ".feature_importance.png"));
            ggplot2::ggsave(imp_plot_file, plot = p_imp,
                            width = 7, height = max(4, top_n * 0.3 + 1),
                            dpi = 300, units = "in");
            msg("    Feature importance plot saved  : ", imp_plot_file);
        }
    }

    msg("\n", strrep("=", 60));
    msg("  Done.");
    msg(strrep("=", 60));

    # ── Return ───────────────────────────────────────────────────────────────
    return(invisible(list(
        tuned_params       = best_params,
        lovo_pred          = lovo_pred,
        metrics            = metrics,
        feature_importance = feature_importance,
        variant_summary    = if (save_plots) variant_summary else NULL
    )));
}

classification_model <- function(
    data,                        # Data frame containing all variables
    predictor,                   # Character vector of predictor variable names
    response,                    # Character string of response variable name (factor or character)
    model        = c("en", "rf", "xgb"),  # Model type
    id_col       = "Record_ID",
    variant_col  = "Variant_name",
    filepath     = NULL,
    model_name   = NULL,
    tune_folds   = 5,
    tune_repeats = 3,
    seed         = 123,
    n_cores      = 1,
    save_plots   = TRUE,
    verbose      = TRUE
    ) {

    # ── Validate arguments ────────────────────────────────────────────────────
    model <- match.arg(model);
    if (is.null(model_name)) model_name <- paste0("lovo_", model, "_clf");
    if (is.null(filepath))   filepath   <- getwd();
    if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE);

    msg <- function(...) if (verbose) message(...);

    # ── Load required packages ────────────────────────────────────────────────
    required_pkgs <- c("caret", "ggplot2", "dplyr");
    model_pkgs <- switch(model,
        en  = c("glmnet"),
        rf  = c("randomForest"),
        xgb = c("xgboost", "Matrix")
    );
    for (pkg in c(required_pkgs, model_pkgs)) {
        if (!requireNamespace(pkg, quietly = TRUE))
            stop("Package '", pkg, "' is required but not installed.");
    }
    library(caret, quietly = TRUE);

    # ── Parallel backend ──────────────────────────────────────────────────────
    if (n_cores > 1) {
        if (!requireNamespace("doParallel", quietly = TRUE))
            stop("Package 'doParallel' is required for parallel computing.");
        cl <- parallel::makeCluster(n_cores);
        doParallel::registerDoParallel(cl);
        on.exit(parallel::stopCluster(cl), add = TRUE);
        msg("  Using ", n_cores, " parallel cores.");
    }

    # ── Prepare data ──────────────────────────────────────────────────────────
    msg("=", strrep("=", 59));
    msg("  Model : ", toupper(model), " (Classification)");
    msg("=", strrep("=", 59));
    msg("\n[1] Preparing data ...");

    if (!variant_col %in% colnames(data))
        stop("Column '", variant_col, "' not found in data.");

    keep_cols <- c(id_col, variant_col, predictor, response);
    temp      <- data[, intersect(keep_cols, colnames(data)), drop = FALSE];
    complete  <- complete.cases(temp[, c(predictor, response)]);
    temp      <- temp[complete, ];
    n         <- nrow(temp);

    # Convert response to factor
    y <- factor(temp[[response]]);
    classes   <- levels(y);
    n_classes <- length(classes);
    is_binary <- n_classes == 2;

    msg("    Samples after NA removal : ", n);
    msg("    Predictors               : ", length(predictor));
    msg("    Classes                  : ", paste(classes, collapse = ", "));
    msg("    Class distribution       : ",
        paste(paste0(classes, "=", table(y)), collapse = ", "));

    unique_variants <- unique(temp[[variant_col]]);
    n_variants      <- length(unique_variants);
    msg("    Unique variants          : ", n_variants);

    variant_sizes <- table(temp[[variant_col]]);
    msg("    Samples per variant      : min=", min(variant_sizes),
        ", median=", median(variant_sizes),
        ", max=", max(variant_sizes));

    sample_ids <- temp[[id_col]];
    X_raw      <- temp[, predictor, drop = FALSE];

    # ── Encode categorical variables ──────────────────────────────────────────
    cat_cols <- predictor[sapply(X_raw, function(v) is.factor(v) || is.character(v))];
    num_cols <- setdiff(predictor, cat_cols);

    if (length(cat_cols) > 0) {
        msg("    Categorical predictors   : ", paste(cat_cols, collapse = ", "));
        dummy_formula <- as.formula(paste("~", paste(cat_cols, collapse = " + ")));
        dummy_mat     <- model.matrix(dummy_formula, data = X_raw)[, -1, drop = FALSE];
        X_proc        <- cbind(X_raw[, num_cols, drop = FALSE], dummy_mat);
    } else {
        X_proc <- X_raw;
    }
    X_proc <- as.data.frame(lapply(X_proc, as.numeric));
    msg("    Final feature count      : ", ncol(X_proc));

    # ── Step 1: Global parameter tuning ──────────────────────────────────────
    msg("\n[2] Global parameter tuning (", tune_folds, "-fold CV x", tune_repeats, " repeats) ...");

    set.seed(seed);
    train_ctrl <- trainControl(
        method          = "repeatedcv",
        number          = tune_folds,
        repeats         = tune_repeats,
        classProbs      = TRUE,
        summaryFunction = if (is_binary) twoClassSummary else multiClassSummary,
        allowParallel   = (n_cores > 1)
    );

    # caret requires syntactically valid factor levels
    y_caret <- factor(make.names(as.character(y)));
    classes_caret <- levels(y_caret);

    if (model == "en") {
        tune_grid <- expand.grid(
            alpha  = seq(0, 1, by = 0.2),
            lambda = 10^seq(-4, 1, length.out = 20)
        );
        set.seed(seed);
        tuned <- caret::train(
            x          = X_proc,
            y          = y_caret,
            method     = "glmnet",
            trControl  = train_ctrl,
            tuneGrid   = tune_grid,
            preProcess = c("center", "scale"),
            metric     = if (is_binary) "ROC" else "Accuracy"
        );

    } else if (model == "rf") {
        tune_grid <- data.frame(
            mtry = unique(round(seq(2, max(2, ncol(X_proc)), length.out = 8)))
        );
        set.seed(seed);
        tuned <- caret::train(
            x         = X_proc,
            y         = y_caret,
            method    = "rf",
            trControl = train_ctrl,
            tuneGrid  = tune_grid,
            ntree     = 500,
            metric    = if (is_binary) "ROC" else "Accuracy"
        );

    } else if (model == "xgb") {
        tune_grid <- expand.grid(
            nrounds          = c(100, 200, 300),
            max_depth        = c(3, 5, 7),
            eta              = c(0.01, 0.05, 0.1),
            gamma            = 0,
            colsample_bytree = c(0.7, 1.0),
            min_child_weight = 1,
            subsample        = 0.8
        );
        set.seed(seed);
        tuned <- caret::train(
            x         = X_proc,
            y         = y_caret,
            method    = "xgbTree",
            trControl = train_ctrl,
            tuneGrid  = tune_grid,
            verbosity = 0,
            metric    = if (is_binary) "ROC" else "Accuracy"
        );
    }

    best_params <- tuned$bestTune;
    msg("    Best parameters:");
    for (nm in names(best_params))
        msg("      ", nm, " = ", best_params[[nm]]);

    if (save_plots) {
        tune_plot_file <- file.path(filepath, paste0(model_name, ".tuning.png"));
        png(tune_plot_file, width = 7, height = 5, units = "in", res = 300);
        print(plot(tuned, main = paste(toupper(model), "- Parameter Tuning")));
        dev.off();
        msg("    Tuning plot saved : ", tune_plot_file);
    }

    # ── Step 2: LOVO-CV ───────────────────────────────────────────────────────
    msg("\n[3] Running Leave-One-Variant-Out CV (", n_variants, " variants, ",
        n, " samples) ...");

    lovo_pred <- data.frame(
        SampleID         = sample_ids,
        Variant          = temp[[variant_col]],
        Actual           = as.character(y),
        Predicted        = NA_character_,
        stringsAsFactors = FALSE
    );

    # Add probability columns
    for (cls in classes) {
        lovo_pred[[paste0("Prob_", cls)]] <- NA_real_;
    }

    # Helper: train and predict
    .fit_predict <- function(train_X, train_y, test_X) {
        train_y_caret <- factor(make.names(as.character(train_y)));

        if (model == "en") {
            library(glmnet, quietly = TRUE);
            pre <- caret::preProcess(train_X, method = c("center", "scale"));
            trX <- predict(pre, train_X);
            teX <- predict(pre, test_X);
            fit <- glmnet::glmnet(
                x      = as.matrix(trX),
                y      = train_y_caret,
                family = if (is_binary) "binomial" else "multinomial",
                alpha  = best_params$alpha,
                lambda = best_params$lambda
            );
            if (is_binary) {
                prob_pos <- as.numeric(predict(fit, newx = as.matrix(teX),
                                               s = best_params$lambda, type = "response"));
                probs <- cbind(1 - prob_pos, prob_pos);
                colnames(probs) <- classes_caret;
            } else {
                probs <- predict(fit, newx = as.matrix(teX),
                                 s = best_params$lambda, type = "response")[, , 1];
            }
            pred_class <- classes_caret[apply(probs, 1, which.max)];

        } else if (model == "rf") {
            library(randomForest, quietly = TRUE);
            fit <- randomForest::randomForest(
                x     = train_X,
                y     = train_y_caret,
                ntree = 500,
                mtry  = best_params$mtry
            );
            probs      <- predict(fit, newdata = test_X, type = "prob");
            pred_class <- as.character(predict(fit, newdata = test_X));

        } else if (model == "xgb") {
            library(xgboost, quietly = TRUE);
            y_int  <- as.integer(train_y_caret) - 1L;
            dtrain <- xgboost::xgb.DMatrix(data = as.matrix(train_X), label = y_int);
            dtest  <- xgboost::xgb.DMatrix(data = as.matrix(test_X));
            fit <- xgboost::xgb.train(
                params = list(
                    objective        = if (is_binary) "binary:logistic" else "multi:softprob",
                    num_class        = if (is_binary) NULL else n_classes,
                    max_depth        = best_params$max_depth,
                    eta              = best_params$eta,
                    gamma            = best_params$gamma,
                    colsample_bytree = best_params$colsample_bytree,
                    min_child_weight = best_params$min_child_weight,
                    subsample        = best_params$subsample,
                    verbosity        = 0
                ),
                data    = dtrain,
                nrounds = best_params$nrounds
            );
            raw <- as.numeric(predict(fit, dtest));
            if (is_binary) {
                probs <- cbind(1 - raw, raw);
                colnames(probs) <- classes_caret;
            } else {
                probs <- matrix(raw, ncol = n_classes, byrow = TRUE);
                colnames(probs) <- classes_caret;
            }
            pred_class <- classes_caret[apply(probs, 1, which.max)];
        }

        # Map back from make.names to original class labels
        class_map  <- setNames(classes, classes_caret);
        pred_orig  <- class_map[pred_class];
        colnames(probs) <- classes;

        return(list(pred = pred_orig, probs = probs));
    }

    # ── LOVO-CV loop ──────────────────────────────────────────────────────────
    for (v_idx in seq_len(n_variants)) {
        variant_name <- unique_variants[v_idx];
        test_mask    <- temp[[variant_col]] == variant_name;
        train_mask   <- !test_mask;
        n_train      <- sum(train_mask);

        if (verbose && (v_idx %% 20 == 0 || v_idx == 1 || v_idx == n_variants))
            msg("    Variant ", v_idx, " / ", n_variants, " : '", variant_name, "'",
                "  (train=", n_train, ", test=", sum(test_mask), ")");

        if (n_train < 5) {
            warning("Variant '", variant_name, "' skipped: only ", n_train,
                    " training samples.");
            next;
        }

        # Check all classes represented in training set
        train_classes <- unique(as.character(y[train_mask]));
        if (length(train_classes) < n_classes) {
            warning("Variant '", variant_name, "' skipped: not all classes in training set.");
            next;
        }

        train_X <- X_proc[train_mask, , drop = FALSE];
        train_y <- y[train_mask];
        test_X  <- X_proc[test_mask,  , drop = FALSE];

        set.seed(seed + v_idx);
        result <- .fit_predict(train_X, train_y, test_X);

        lovo_pred$Predicted[test_mask] <- result$pred;
        for (cls in classes) {
            lovo_pred[[paste0("Prob_", cls)]][test_mask] <- result$probs[, cls];
        }
    }

    # ── Step 3: Performance metrics ───────────────────────────────────────────
    msg("\n[4] Computing performance metrics ...");

    valid_idx  <- !is.na(lovo_pred$Predicted);
    n_valid    <- sum(valid_idx);
    msg("    Samples with predictions : ", n_valid, " / ", n);

    actual_v    <- factor(lovo_pred$Actual[valid_idx],    levels = classes);
    predicted_v <- factor(lovo_pred$Predicted[valid_idx], levels = classes);

    cm       <- caret::confusionMatrix(predicted_v, actual_v);
    accuracy <- as.numeric(cm$overall["Accuracy"]);
    kappa    <- as.numeric(cm$overall["Kappa"]);

    # AUC (binary only)
    auc_val <- NA_real_;
    if (is_binary && requireNamespace("pROC", quietly = TRUE)) {
        prob_col <- paste0("Prob_", classes[2]);
        auc_val  <- as.numeric(pROC::auc(
            pROC::roc(actual_v, lovo_pred[[prob_col]][valid_idx], quiet = TRUE)
        ));
    }

    metrics <- data.frame(
        Model    = model,
        CV_type  = "LOVO",
        N        = n_valid,
        Variants = n_variants,
        Accuracy = round(accuracy, 4),
        Kappa    = round(kappa,    4),
        AUC      = round(auc_val,  4)
    );

    msg("    Accuracy = ", round(accuracy, 4));
    msg("    Kappa    = ", round(kappa,    4));
    if (!is.na(auc_val)) msg("    AUC      = ", round(auc_val, 4));
    msg("\n    Confusion Matrix:");
    print(cm$table);

    # ── Step 4: Feature importance on full dataset ────────────────────────────
    msg("\n[5] Computing feature importance on full dataset ...");
    feature_importance <- NULL;
    y_caret_full <- factor(make.names(as.character(y)));

    if (model == "en") {
        library(glmnet, quietly = TRUE);
        pre_full <- caret::preProcess(X_proc, method = c("center", "scale"));
        X_scaled <- predict(pre_full, X_proc);
        fit_full <- glmnet::glmnet(
            x      = as.matrix(X_scaled),
            y      = y_caret_full,
            family = if (is_binary) "binomial" else "multinomial",
            alpha  = best_params$alpha,
            lambda = best_params$lambda
        );
        coefs <- coef(fit_full, s = best_params$lambda);
        if (is_binary) {
            coef_mat <- as.matrix(coefs)[-1, , drop = FALSE];
            feature_importance <- data.frame(
                Feature     = rownames(coef_mat),
                Importance  = abs(coef_mat[, 1]),
                Coefficient = coef_mat[, 1]
            );
        } else {
            coef_mat <- do.call(cbind, lapply(coefs, function(x) as.matrix(x)[-1, ]));
            colnames(coef_mat) <- classes;
            feature_importance <- data.frame(
                Feature    = rownames(coef_mat),
                Importance = rowMeans(abs(coef_mat))
            );
            feature_importance <- cbind(feature_importance, as.data.frame(coef_mat));
        }
        feature_importance <- feature_importance[
            order(feature_importance$Importance, decreasing = TRUE), ];

    } else if (model == "rf") {
        library(randomForest, quietly = TRUE);
        fit_full <- randomForest::randomForest(
            x          = X_proc,
            y          = y_caret_full,
            ntree      = 500,
            mtry       = best_params$mtry,
            importance = TRUE
        );
        imp_mat <- randomForest::importance(fit_full);
        feature_importance <- data.frame(
            Feature    = rownames(imp_mat),
            Importance = imp_mat[, "MeanDecreaseAccuracy"],
            MeanDecreaseGini = imp_mat[, "MeanDecreaseGini"]
        );
        feature_importance <- feature_importance[
            order(feature_importance$Importance, decreasing = TRUE), ];

    } else if (model == "xgb") {
        library(xgboost, quietly = TRUE);
        y_int_full  <- as.integer(y_caret_full) - 1L;
        dtrain_full <- xgboost::xgb.DMatrix(data = as.matrix(X_proc), label = y_int_full);
        fit_full <- xgboost::xgb.train(
            params = list(
                objective        = if (is_binary) "binary:logistic" else "multi:softprob",
                num_class        = if (is_binary) NULL else n_classes,
                max_depth        = best_params$max_depth,
                eta              = best_params$eta,
                gamma            = best_params$gamma,
                colsample_bytree = best_params$colsample_bytree,
                min_child_weight = best_params$min_child_weight,
                subsample        = best_params$subsample,
                verbosity        = 0
            ),
            data    = dtrain_full,
            nrounds = best_params$nrounds
        );
        imp_raw <- xgboost::xgb.importance(feature_names = colnames(X_proc), model = fit_full);
        feature_importance <- as.data.frame(imp_raw);
        colnames(feature_importance)[1] <- "Feature";
        colnames(feature_importance)[2] <- "Importance";
    }

    # ── Step 5: Save outputs ──────────────────────────────────────────────────
    msg("\n[6] Saving outputs to : ", filepath);

    pred_file <- file.path(filepath, paste0(model_name, ".lovo_predictions.txt"));
    write.table(lovo_pred, file = pred_file, row.names = FALSE, quote = FALSE, sep = "\t");
    msg("    Predictions saved  : ", pred_file);

    metric_file <- file.path(filepath, paste0(model_name, ".metrics.txt"));
    write.table(metrics, file = metric_file, row.names = FALSE, quote = FALSE, sep = "\t");
    msg("    Metrics saved      : ", metric_file);

    if (!is.null(feature_importance)) {
        imp_file <- file.path(filepath, paste0(model_name, ".feature_importance.txt"));
        write.table(feature_importance, file = imp_file, row.names = FALSE, quote = FALSE, sep = "\t");
        msg("    Importance saved   : ", imp_file);
    }

    if (save_plots) {
        # (i) Confusion matrix heatmap
        cm_df <- as.data.frame(cm$table);
        colnames(cm_df) <- c("Predicted", "Actual", "Freq");
        cm_file <- file.path(filepath, paste0(model_name, ".confusion_matrix.png"));
        p_cm <- ggplot2::ggplot(cm_df, ggplot2::aes(x = Actual, y = Predicted, fill = Freq)) +
            ggplot2::geom_tile() +
            ggplot2::geom_text(ggplot2::aes(label = Freq), size = 5) +
            ggplot2::scale_fill_gradient(low = "white", high = "#2166AC") +
            ggplot2::labs(title = paste(toupper(model), "- Confusion Matrix"),
                          x = "Actual", y = "Predicted") +
            ggplot2::theme_bw(base_size = 12);
        ggplot2::ggsave(cm_file, p_cm, width = 5, height = 4, dpi = 300);
        msg("    Confusion matrix plot saved : ", cm_file);

        # (ii) ROC curve (binary only)
        if (is_binary && requireNamespace("pROC", quietly = TRUE)) {
            prob_col <- paste0("Prob_", classes[2]);
            roc_obj  <- pROC::roc(actual_v, lovo_pred[[prob_col]][valid_idx], quiet = TRUE);
            roc_df   <- data.frame(
                FPR = 1 - roc_obj$specificities,
                TPR = roc_obj$sensitivities
            );
            roc_file <- file.path(filepath, paste0(model_name, ".roc_curve.png"));
            p_roc <- ggplot2::ggplot(roc_df, ggplot2::aes(x = FPR, y = TPR)) +
                ggplot2::geom_line(color = "#2166AC", linewidth = 1) +
                ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
                ggplot2::annotate("text", x = 0.75, y = 0.25,
                                  label = paste0("AUC = ", round(auc_val, 3)),
                                  size = 5, fontface = "bold") +
                ggplot2::labs(title = paste(toupper(model), "- ROC Curve"),
                              x = "False Positive Rate", y = "True Positive Rate") +
                ggplot2::theme_bw(base_size = 12);
            ggplot2::ggsave(roc_file, p_roc, width = 5, height = 5, dpi = 300);
            msg("    ROC curve saved             : ", roc_file);
        }

        # (iii) Feature importance bar plot (top 20)
        if (!is.null(feature_importance) && nrow(feature_importance) > 0) {
            top_n  <- min(20, nrow(feature_importance));
            top_df <- feature_importance[seq_len(top_n), ];
            top_df$Feature <- factor(top_df$Feature, levels = rev(top_df$Feature));
            p_imp <- ggplot2::ggplot(top_df, ggplot2::aes(x = Feature, y = Importance)) +
                ggplot2::geom_bar(stat = "identity", fill = "#4393C3") +
                ggplot2::coord_flip() +
                ggplot2::labs(title = paste(toupper(model), "- Feature Importance (Top", top_n, ")"),
                              x = NULL, y = "Importance") +
                ggplot2::theme_bw(base_size = 11);
            imp_plot_file <- file.path(filepath, paste0(model_name, ".feature_importance.png"));
            ggplot2::ggsave(imp_plot_file, p_imp,
                            width = 7, height = max(4, top_n * 0.3 + 1), dpi = 300);
            msg("    Feature importance plot saved : ", imp_plot_file);
        }
    }

    msg("\n", strrep("=", 60));
    msg("  Done.");
    msg(strrep("=", 60));

    return(invisible(list(
        tuned_params       = best_params,
        lovo_pred          = lovo_pred,
        metrics            = metrics,
        confusion_matrix   = cm,
        feature_importance = feature_importance
    )));
}

model_evaluation <- function(
        model,
        data,
        type            = c("regression", "classification"),
        age_name        = NULL,
        age_cutoff      = NULL,
        residual_cutoff = 10,
        axis_label      = NULL,
        filepath        = NULL
    ) {

    type <- match.arg(type);
    if (is.null(filepath)) filepath <- getwd();
    if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE);

    select_colnames <- c("Record_ID", "Age_VABS", "Age_VABS_1st", "Age_VABS_Diff",
                         "VABS_ABC", "VABS_ABC_Diff", "VABS_ABC_Rate",
                         "Variant_name", "GBTM_cluster");
    select_colnames <- intersect(select_colnames, colnames(data));

    model_names <- names(model);

    if (type == "regression") {

        # ── 1a. Actual vs Predicted plots (regression) ────────────────────────
        for (model_name in model_names) {
            m <- model[[model_name]];

            temp_data <- m$loocv_pred;
            temp_data$Residual <- abs(temp_data$Actual - temp_data$Predicted);

            temp_data <- merge(temp_data, data[, select_colnames],
                               by.x = "SampleID", by.y = "Record_ID",
                               all.x = TRUE, all.y = FALSE);

            temp_data[temp_data$Residual <= residual_cutoff, ]$Variant_name <- "";

            if (!is.null(age_name) & !is.null(age_cutoff)) {
                temp_data <- temp_data[!is.na(temp_data[[age_name]]) &
                                       temp_data[[age_name]] >= age_cutoff, ];
            };

            x_range  <- range(c(temp_data$Actual, temp_data$Predicted), na.rm = TRUE);
            y_range  <- x_range;
            r2_value <- cor(temp_data$Actual, temp_data$Predicted, use = "complete.obs") ^ 2;

            p <- ggplot2::ggplot(temp_data, ggplot2::aes(x = Actual, y = Predicted)) +
                ggplot2::theme_classic() +
                ggplot2::xlim(x_range) + ggplot2::ylim(y_range) +
                ggplot2::geom_point(
                    ggplot2::aes(
                        fill  = Age_VABS,
                        color = ifelse(is.na(GBTM_cluster), NA, as.character(GBTM_cluster))
                    ),
                    shape = 21, size = 4, alpha = 0.8, stroke = 1.2
                ) +
                ggplot2::scale_fill_gradient(low = "#98DF8A", high = "#FF9896",
                                             name = "Age at endpoint") +
                ggplot2::scale_color_manual(
                    values   = c("1" = "#17BECF", "2" = "#FF7F0E"),
                    na.value = "grey60",
                    name     = "Patient cluster",
                    breaks   = c("1", "2")
                ) +
                ggplot2::guides(
                    color = ggplot2::guide_legend(
                        override.aes = list(shape = 21, fill = "grey70", size = 4)
                    ),
                    fill = ggplot2::guide_colorbar()
                ) +
                ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
                ggrepel::geom_text_repel(
                    ggplot2::aes(label = Variant_name),
                    size         = 3,
                    box.padding  = 0.5,
                    max.overlaps = 30
                ) +
                ggplot2::annotate("text",
                                  x     = -Inf, y = Inf,
                                  label = paste("Pearson~r^2 == ", sprintf("%0.2f", r2_value), sep = ""),
                                  parse = TRUE,
                                  size  = 5, fontface = "bold",
                                  hjust = -0.2, vjust = 1.5) +
                ggplot2::labs(
                    title = NULL,
                    x     = paste("Actual", axis_label),
                    y     = paste("Predicted", axis_label)
                ) +
                ggplot2::theme(legend.position = c(0.90, 0.2));

            filename <- file.path(filepath, paste0("Actual_vs_Predicted.", model_name, ".png"));
            ggplot2::ggsave(filename, p, width = 10, height = 8, units = "in", dpi = 300);
            message("Saved: ", filename);
        };

        # ── 2a. Ablation analysis (regression) ────────────────────────────────
        metrics_list <- lapply(model_names, function(model_name) {
            df <- model[[model_name]]$metrics;
            df$Model_config <- model_name;
            return(df);
        });

        ablation_table <- do.call(rbind, metrics_list);
        ablation_table <- ablation_table[, c("Model_config",
                                             setdiff(colnames(ablation_table), "Model_config"))];

        if ("Full" %in% model_names) {
            full_metrics <- ablation_table[ablation_table$Model_config == "Full", ];
            for (metric in c("RMSE", "MAE", "R2", "Pearson")) {
                if (metric %in% colnames(ablation_table)) {
                    ablation_table[[paste0("Delta_", metric)]] <- round(
                        ablation_table[[metric]] - full_metrics[[metric]], 4
                    );
                };
            };
        };

        ablation_file <- file.path(filepath, "Ablation_analysis.txt");
        write.table(ablation_table, ablation_file,
                    sep = "\t", quote = FALSE, row.names = FALSE);
        message("Ablation table saved: ", ablation_file);

    } else {

        # ── 1b. Confusion matrix plots (classification) ───────────────────────
        for (model_name in model_names) {
            m  <- model[[model_name]];
            cm <- m$confusion_matrix;

            cm_df <- as.data.frame(cm$table);
            colnames(cm_df) <- c("Predicted", "Actual", "Freq");

            # Annotate with row percentage
            cm_df <- cm_df %>%
                dplyr::group_by(Actual) %>%
                dplyr::mutate(Pct = round(100 * Freq / sum(Freq), 1)) %>%
                dplyr::ungroup() %>%
                dplyr::mutate(label = paste0(Freq, "\n(", Pct, "%)"));

            p <- ggplot2::ggplot(cm_df, ggplot2::aes(x = Actual, y = Predicted, fill = Freq)) +
                ggplot2::geom_tile(color = "white") +
                ggplot2::geom_text(ggplot2::aes(label = label), size = 4.5) +
                ggplot2::scale_fill_gradient(low = "white", high = "#2166AC",
                                             name = "Count") +
                ggplot2::labs(title = NULL, x = "Actual", y = "Predicted") +
                ggplot2::theme_classic(base_size = 13);

            filename <- file.path(filepath, paste0("Confusion_matrix.", model_name, ".png"));
            ggplot2::ggsave(filename, p, width = 6, height = 5, units = "in", dpi = 300);
            message("Saved: ", filename);

            # ROC curve if binary and prob columns available
            prob_cols <- grep("^Prob_", colnames(m$lovo_pred), value = TRUE);
            if (length(prob_cols) == 2 && requireNamespace("pROC", quietly = TRUE)) {
                valid_idx <- !is.na(m$lovo_pred$Predicted);
                actual_v  <- factor(m$lovo_pred$Actual[valid_idx]);
                prob_pos  <- m$lovo_pred[[prob_cols[2]]][valid_idx];
                roc_obj   <- pROC::roc(actual_v, prob_pos, quiet = TRUE);
                auc_val   <- as.numeric(pROC::auc(roc_obj));
                roc_df    <- data.frame(
                    FPR = 1 - roc_obj$specificities,
                    TPR = roc_obj$sensitivities
                );

                p_roc <- ggplot2::ggplot(roc_df, ggplot2::aes(x = FPR, y = TPR)) +
                    ggplot2::geom_line(color = "#2166AC", linewidth = 1) +
                    ggplot2::geom_abline(intercept = 0, slope = 1,
                                         linetype = "dashed", color = "gray50") +
                    ggplot2::annotate("text",
                                      x = -Inf, y = -Inf,
                                      label = paste0("AUC = ", round(auc_val, 3)),
                                      size = 5, fontface = "bold",
                                      hjust = -0.2, vjust = -1) +
                    ggplot2::labs(title = NULL,
                                  x = "False Positive Rate",
                                  y = "True Positive Rate") +
                    ggplot2::theme_classic(base_size = 13);

                roc_file <- file.path(filepath, paste0("ROC_curve.", model_name, ".png"));
                ggplot2::ggsave(roc_file, p_roc, width = 5, height = 5, units = "in", dpi = 300);
                message("Saved: ", roc_file);
            };
        };

        # ── 2b. Ablation analysis (classification) ────────────────────────────
        metrics_list <- lapply(model_names, function(model_name) {
            df <- model[[model_name]]$metrics;
            df$Model_config <- model_name;
            return(df);
        });

        ablation_table <- do.call(rbind, metrics_list);
        ablation_table <- ablation_table[, c("Model_config",
                                             setdiff(colnames(ablation_table), "Model_config"))];

        if ("Full" %in% model_names) {
            full_metrics <- ablation_table[ablation_table$Model_config == "Full", ];
            for (metric in c("Accuracy", "Kappa", "AUC")) {
                if (metric %in% colnames(ablation_table)) {
                    ablation_table[[paste0("Delta_", metric)]] <- round(
                        ablation_table[[metric]] - full_metrics[[metric]], 4
                    );
                };
            };
        };

        ablation_file <- file.path(filepath, "Ablation_analysis.txt");
        write.table(ablation_table, ablation_file,
                    sep = "\t", quote = FALSE, row.names = FALSE);
        message("Ablation table saved: ", ablation_file);
    };
};
# === Function End ===

# =============== #
# Load KIF1A data #
# =============== #

# ------------------------- #
# KIF1A variable annotation #
# ------------------------- #

filename <- paste(project_path, "/KIF1A_anno.xlsx", sep="");
KIF1A_anno <- read.xlsx(filename, sheet="KIF1A_anno.v3", sep.names=" ");

# ------------------- #
# KIF1A variants data #
# ------------------- #

filename <- paste(rawdata_path, "/clinical_data/Deidentified KIF1A Data November 2025.xlsx", sep="");
KIF1A_variant <- read.xlsx(filename, sheet = "Cohort Overview", startRow = 2);
KIF1A_variant[is.na(KIF1A_variant$Record.ID), ]$Record.ID <- "K_024";

colnames(KIF1A_variant) <- c(
    "Record_ID", "K_ID", "Birth_date", "Age", "Sex", "Country", "Race", "Ethnicity", "Mutation_type", "c_dot", 
    "p_dot", "Transcript", "Remove_1", "Remove_2", "Remove_3", "Remove_4", "Inheritance", "ACMG_classification", 
    "Bugrahan_classification", "Historical_MHI_Y1", "Historical_MHI_Y2", "Historical_MHI_Y3", "Current_MHI_Y1", 
    "Current_MHI_Y2", "ASCEND", "KOALA_Y1", "KOALA_Y2", "KOALA_Y3", "KOALA_Y4"
);

KIF1A_variant <- KIF1A_variant[, colnames(KIF1A_variant)[!grepl("^Remove", colnames(KIF1A_variant))]];

KIF1A_variant$Death <- "No";
KIF1A_variant[grepl("^DECEASED", KIF1A_variant$Record_ID), ]$Death <- "Yes";
KIF1A_variant$Age <- round(as.numeric(KIF1A_variant$Age), 2);

# === Correct variant data ===
KIF1A_variant$c_dot <- trimws(KIF1A_variant$c_dot);
for (i in c(1:nrow(KIF1A_variant))) {
    if (grepl(">", KIF1A_variant$c_dot[i])) {
        KIF1A_variant$c_dot[i] <- gsub(" ", "", KIF1A_variant$c_dot[i]);
    }
}

# Correct transcript reference
KIF1A_variant$Transcript <- "NM_001244008.2";

# For cases carry multiple variants, keep the most severe one
KIF1A_variant[KIF1A_variant$Record_ID == "NHS01-225", "c_dot"] <- "c.536A>G";
KIF1A_variant[KIF1A_variant$Record_ID == "NHS01-225", "p_dot"] <- "p.Glu179Gly";
KIF1A_variant[KIF1A_variant$Record_ID == "NHS01-229", "c_dot"] <- "c.685C>T";
KIF1A_variant[KIF1A_variant$Record_ID == "NHS01-229", "p_dot"] <- "p.Arg229Cys";
KIF1A_variant[KIF1A_variant$Record_ID == "NHS01-285", "c_dot"] <- "c.2664del";
KIF1A_variant[KIF1A_variant$Record_ID == "NHS01-285", "p_dot"] <- "p.Thr889ProfsTer62";

# Correct c_dot
KIF1A_variant[KIF1A_variant$Record_ID == "NHS01-170", "c_dot"] <- "c.835G>C";
KIF1A_variant[KIF1A_variant$Record_ID == "NHS01-244", "c_dot"] <- "c.1451C>A";
KIF1A_variant[KIF1A_variant$Record_ID == "NHS01-215", "c_dot"] <- "c.4711A>T";

# Create de novo data for annotation
KIF1A_variant$ID <- paste(KIF1A_variant$Transcript, KIF1A_variant$c_dot, sep = ":");
KIF1A_variant[!grepl("^c.", KIF1A_variant$c_dot), "ID"] <- "";

unique_ID <- unique(KIF1A_variant[KIF1A_variant$ID != "", ]$ID);
filename <- paste(rawdata_path, "/predicted_score/KIF1A_variants_for_VEP_anno.txt", sep="");
writeLines(unique_ID, filename);

# ------------------- #
# KIF1A clinical data #
# ------------------- #

# === VABS score ===
filename <- paste(rawdata_path, "/clinical_data/Deidentified KIF1A Data November 2025.xlsx", sep="");
KIF1A_VABS <- read.xlsx(filename, sheet = "Vinelands", startRow = 1);
colnames(KIF1A_VABS) <- c("Record_ID", "Instance", "Age_VABS", "VABS_ABC", "VABS_communication", "VABS_DLS", "VABS_motor", "VABS_socialization");
KIF1A_VABS$Key <- paste(KIF1A_VABS$Record_ID, KIF1A_VABS$Instance, sep = "-");

# For cases with NULL motor score, assign the bottom-truncated value
KIF1A_VABS[is.na(KIF1A_VABS$VABS_motor), "VABS_motor"] <- 20;

# === Clinical data ===
filename <- paste(rawdata_path, "/clinical_data/Deidentified KIF1A Data November 2025.xlsx", sep="");
KIF1A_clinical <- read.xlsx(filename, sheet = "Yufeng variables", startRow = 3);

# Remove the blank column "X42" in the clinical data
KIF1A_clinical <- KIF1A_clinical[, setdiff(colnames(KIF1A_clinical), "X42")];

# Rename the columns
clinical_colnames <- c(
    "Age_survey", "Date_survey", "Age_first_concern", "Age_measured", "Height", "Weight", "Vision_concern", "Optic_nerve_change", 
    "Age_optic_nerve_change", "Seizure_any", "Age_seizure_onset", "Seizure_frequency_current", "Seizure_medication", "Seizure_resistant", 
    "Imaging_any", "CT_done", "MRI_done", "CT_result", "CT_result_description", "MRI_result", "MRI_result_description", "EEG_done", 
    "EEG_result", "EEG_result_description", "Neurologic_diagnosis", "Hypertonia", "Age_hypertonia", "Hypertonia_resolved", 
    "Age_hypertonia_resolved", "Hypotonia", "Age_hypotonia", "Hypotonia_resolved", "Age_hypotonia_resolved", "Microcephaly", 
    "Age_microcephaly", "Movement_abnormal", "Ataxia", "Age_ataxia", "Ataxia_resolved", "Dystonia", "Age_dystonia", "Dystonia_resolved", 
    "Neuropathy", "Age_neuropathy", "Neuropathy_resolved", "Developmental_concern", "Regression_hand", "Age_regression_hand", 
    "Trigger_regression_hand", "Duration_regression_hand", "Regression_language", "Age_regression_language", "Trigger_regression_language", 
    "Duration_regression_language", "Regression_daily_living", "Age_regression_daily_living", "Trigger_regression_daily_living", 
    "Duration_regression_daily_living", "Regression_motor", "Age_regression_motor", "Trigger_regression_motor", "Duration_regression_motor", 
    "Age_motor_recovered", "Independent_sit", "Age_sit", "Independent_walk", "Age_walk", "Speech", "Age_first_word"
);
colnames(KIF1A_clinical) <- c("Record_ID", "Instance", clinical_colnames);

# The missing Record_ID in the clinical data is "K_024", which can be inferred from the variants data.
KIF1A_clinical[is.na(KIF1A_clinical$Record_ID), ]$Record_ID <- "K_024";

# Correct the Instance column for each patient, Hx - Historical
KIF1A_clinical$Instance_old <- KIF1A_clinical$Instance;
KIF1A_clinical$Is_Hx <- grepl("^Hx", KIF1A_clinical$Instance_old);
KIF1A_clinical$Num_part <- as.numeric(gsub("^Hx", "", KIF1A_clinical$Instance_old));
KIF1A_clinical <- KIF1A_clinical[order(KIF1A_clinical$Record_ID, -KIF1A_clinical$Is_Hx, KIF1A_clinical$Num_part), ];
KIF1A_clinical$Instance <- ave(KIF1A_clinical$Instance, KIF1A_clinical$Record_ID, FUN = seq_along);
KIF1A_clinical$Instance <- as.numeric(KIF1A_clinical$Instance);
KIF1A_clinical <- KIF1A_clinical[, setdiff(colnames(KIF1A_clinical), c("Instance_old", "Is_Hx", "Num_part"))];

# Merge the clinical data with the VABS data
KIF1A_clinical$Key <- paste(KIF1A_clinical$Record_ID, KIF1A_clinical$Instance, sep = "-");
KIF1A_clinical <- merge(KIF1A_clinical, KIF1A_VABS[, setdiff(colnames(KIF1A_VABS), c("Record_ID", "Instance"))], 
                        by = "Key", all.x = TRUE, all.y = FALSE);
KIF1A_clinical <- KIF1A_clinical[, setdiff(colnames(KIF1A_clinical), "Key")];

# Remove the records with all missing values in the clinical data (except for the baseline information)
# select_colnames <- c("Record_ID", "Instance", "Age_survey", "Date_survey", "Age_first_concern", "Age_measured", "Height", "Weight");
# KIF1A_clinical$NA_ratio <- rowSums(is.na(KIF1A_clinical[, setdiff(colnames(KIF1A_clinical), select_colnames)])) / 
#                            (ncol(KIF1A_clinical) - length(select_colnames));
# KIF1A_clinical <- KIF1A_clinical[KIF1A_clinical$NA_ratio < 1, ];
# KIF1A_clinical <- KIF1A_clinical[, setdiff(colnames(KIF1A_clinical), "NA_ratio")];

# Standardize the clinical data
for (i in c(1:ncol(KIF1A_clinical))) {
    if (any(grepl("^abnormal$", tolower(KIF1A_clinical[, i])))) {
        KIF1A_clinical[grepl("^abnormal$", tolower(KIF1A_clinical[, i])), i] <- "Abnormal";
    }
    if (any(grepl("^normal$", tolower(KIF1A_clinical[, i])))) {
        KIF1A_clinical[grepl("^normal$", tolower(KIF1A_clinical[, i])), i] <- "Normal";
    }
    if (any(grepl("^yes$", tolower(KIF1A_clinical[, i])))) {
        KIF1A_clinical[grepl("^yes$", tolower(KIF1A_clinical[, i])), i] <- "Yes";
    }
    if (any(grepl("^no$", tolower(KIF1A_clinical[, i])))) {
        KIF1A_clinical[grepl("^no$", tolower(KIF1A_clinical[, i])), i] <- "No";
    }
    if (any(grepl("^checked$", tolower(KIF1A_clinical[, i])))) {
        KIF1A_clinical[grepl("^checked$", tolower(KIF1A_clinical[, i])), i] <- "Yes";
    }
    if (any(grepl("^unchecked$", tolower(KIF1A_clinical[, i])))) {
        KIF1A_clinical[grepl("^unchecked$", tolower(KIF1A_clinical[, i])), i] <- "No";
    }
    if (any(grepl("^unknown$", tolower(KIF1A_clinical[, i])))) {
        KIF1A_clinical[grepl("^unknown$", tolower(KIF1A_clinical[, i])), i] <- "Unknown";
    }
    if (any(grepl("^not sure$", tolower(KIF1A_clinical[, i])))) {
        KIF1A_clinical[grepl("^not sure$", tolower(KIF1A_clinical[, i])), i] <- "Unknown";
    }
}

# Correct clinical records for vision
# - If there is any record of the age of optic nerve change, then the patient must have optic nerve change and vision concern
KIF1A_clinical[which(!is.na(KIF1A_clinical$Age_optic_nerve_change)), "Optic_nerve_change"] <- "Yes";
KIF1A_clinical[which(!is.na(KIF1A_clinical$Optic_nerve_change) & 
                     KIF1A_clinical$Optic_nerve_change == "Yes"), "Vision_concern"] <- "Yes";
KIF1A_clinical[which(is.na(KIF1A_clinical$Optic_nerve_change)), "Optic_nerve_change"] <- "Unknown";
KIF1A_clinical[which(is.na(KIF1A_clinical$Vision_concern)), "Vision_concern"] <- "Unknown";

# Correct clinical records for seizure
# - If there is any record of seizure onset age, current seizure frequency, seizure medication, 
#   or drug-resistant seizure, then the patient must have seizure
KIF1A_clinical[which(!is.na(KIF1A_clinical$Age_seizure_onset) | 
                     !is.na(KIF1A_clinical$Seizure_frequency_current) | 
                     !is.na(KIF1A_clinical$Seizure_medication) | 
                     !is.na(KIF1A_clinical$Seizure_resistant)), "Seizure_any"] <- "Yes";
KIF1A_clinical[which(is.na(KIF1A_clinical$Seizure_any) | 
                     KIF1A_clinical$Seizure_any == "Unknown"), "Seizure_any"] <- "Unknown";

# Correct clinical records for brain imaging
KIF1A_clinical[which(!is.na(KIF1A_clinical$CT_result) | !is.na(KIF1A_clinical$CT_result_description)), "CT_done"] <- "Yes";
KIF1A_clinical[which(!is.na(KIF1A_clinical$MRI_result) | !is.na(KIF1A_clinical$MRI_result_description)), "MRI_done"] <- "Yes";
KIF1A_clinical[which((!is.na(KIF1A_clinical$CT_done) & KIF1A_clinical$CT_done == "Yes") | 
                     (!is.na(KIF1A_clinical$MRI_done) & KIF1A_clinical$MRI_done == "Yes")), "Brain_imaging_any"] <- "Yes";
KIF1A_clinical[which((!is.na(KIF1A_clinical$CT_done) & KIF1A_clinical$CT_done == "No") & 
                     (!is.na(KIF1A_clinical$MRI_done) & KIF1A_clinical$MRI_done == "No")), "Brain_imaging_any"] <- "No";
KIF1A_clinical[which(is.na(KIF1A_clinical$Brain_imaging_any)), "Brain_imaging_any"] <- "Unknown";
KIF1A_clinical[which(is.na(KIF1A_clinical$CT_done)), "CT_done"] <- "Unknown";
KIF1A_clinical[which(is.na(KIF1A_clinical$MRI_done)), "MRI_done"] <- "Unknown";

# Correct clinical records for EEG
# - If there is any record of EEG result description, then the patient must have abnormal EEG
KIF1A_clinical[which(!is.na(KIF1A_clinical$EEG_result_description)), "EEG_result"] <- "Abnormal";
KIF1A_clinical[which(!is.na(KIF1A_clinical$EEG_result)), "EEG_done"] <- "Yes";
KIF1A_clinical[which(is.na(KIF1A_clinical$EEG_done)), "EEG_done"] <- "Unknown";
KIF1A_clinical[which(is.na(KIF1A_clinical$EEG_result)), "EEG_result"] <- "Unknown";

# Correct clinical records for hypertonia
# - If there is any record of the age of hypertonia, or whether hypertonia is resolved, 
#   then the patient must have hypertonia
KIF1A_clinical[which(!is.na(KIF1A_clinical$Age_hypertonia) | 
                     !is.na(KIF1A_clinical$Hypertonia_resolved) | 
                     !is.na(KIF1A_clinical$Age_hypertonia_resolved)), "Hypertonia"] <- "Yes";
KIF1A_clinical[which(is.na(KIF1A_clinical$Hypertonia)), "Hypertonia"] <- "Unknown";

# Correct clinical records for hypotonia
# - If there is any record of the age of hypotonia, or whether hypotonia is resolved, 
#   then the patient must have hypotonia
KIF1A_clinical[which(!is.na(KIF1A_clinical$Age_hypotonia) | 
                     !is.na(KIF1A_clinical$Hypotonia_resolved) | 
                     !is.na(KIF1A_clinical$Age_hypotonia_resolved)), "Hypotonia"] <- "Yes";
KIF1A_clinical[which(is.na(KIF1A_clinical$Hypotonia)), "Hypotonia"] <- "Unknown";

# Correct clinical records for microcephaly
KIF1A_clinical[which(!is.na(KIF1A_clinical$Age_microcephaly)), "Microcephaly"] <- "Yes";
KIF1A_clinical[which(is.na(KIF1A_clinical$Microcephaly)), "Microcephaly"] <- "Unknown";

# Correct clinical records for ataxia
KIF1A_clinical[which(!is.na(KIF1A_clinical$Age_ataxia) | 
                     !is.na(KIF1A_clinical$Ataxia_resolved)), "Ataxia"] <- "Yes";
KIF1A_clinical[which(is.na(KIF1A_clinical$Ataxia)), "Ataxia"] <- "Unknown";

# Correct clinical records for dystonia
KIF1A_clinical[which(!is.na(KIF1A_clinical$Age_dystonia) | 
                     !is.na(KIF1A_clinical$Dystonia_resolved)), "Dystonia"] <- "Yes";
KIF1A_clinical[which(is.na(KIF1A_clinical$Dystonia)), "Dystonia"] <- "Unknown";

# Correct clinical records for neuropathy
KIF1A_clinical[which(!is.na(KIF1A_clinical$Age_neuropathy) | 
                     !is.na(KIF1A_clinical$Neuropathy_resolved)), "Neuropathy"] <- "Yes";
KIF1A_clinical[which(is.na(KIF1A_clinical$Neuropathy)), "Neuropathy"] <- "Unknown";

# Correct clinical records for movement abnormalities
# - If there is any record of ataxia, dystonia, or neuropathy, then the patient 
#   must have movement abnormalities
KIF1A_clinical[which(KIF1A_clinical$Ataxia == "Yes" | 
                     KIF1A_clinical$Dystonia == "Yes" | 
                     KIF1A_clinical$Neuropathy == "Yes"), "Movement_abnormal"] <- "Yes";
KIF1A_clinical[which(is.na(KIF1A_clinical$Movement_abnormal)), "Movement_abnormal"] <- "Unknown";

# Correct clinical records for regression of hand skills
KIF1A_clinical[which(!is.na(KIF1A_clinical$Age_regression_hand) | 
                     !is.na(KIF1A_clinical$Trigger_regression_hand) | 
                     !is.na(KIF1A_clinical$Duration_regression_hand)), "Regression_hand"] <- "Yes";
KIF1A_clinical[which(is.na(KIF1A_clinical$Regression_hand)), "Regression_hand"] <- "Unknown";

# Correct clinical records for regression of language
KIF1A_clinical[which(!is.na(KIF1A_clinical$Age_regression_language) | 
                     !is.na(KIF1A_clinical$Trigger_regression_language) | 
                     !is.na(KIF1A_clinical$Duration_regression_language)), "Regression_language"] <- "Yes";
KIF1A_clinical[which(is.na(KIF1A_clinical$Regression_language)), "Regression_language"] <- "Unknown";

# Correct clinical records for regression of daily living
KIF1A_clinical[which(!is.na(KIF1A_clinical$Age_regression_daily_living) | 
                     !is.na(KIF1A_clinical$Trigger_regression_daily_living) | 
                     !is.na(KIF1A_clinical$Duration_regression_daily_living)), "Regression_daily_living"] <- "Yes";
KIF1A_clinical[which(is.na(KIF1A_clinical$Regression_daily_living)), "Regression_daily_living"] <- "Unknown";

# Correct clinical records for regression gross motor
KIF1A_clinical[which(!is.na(KIF1A_clinical$Age_regression_motor) | 
                     !is.na(KIF1A_clinical$Trigger_regression_motor) | 
                     !is.na(KIF1A_clinical$Duration_regression_motor) | 
                     !is.na(KIF1A_clinical$Age_motor_recovered)), "Regression_motor"] <- "Yes";
KIF1A_clinical[which(is.na(KIF1A_clinical$Regression_motor)), "Regression_motor"] <- "Unknown";

# Correct clinical records for developmental concerns
# - If there is any record of regression in hand skills, language, daily living, or motor skills,
#   then the patient must have developmental concerns
KIF1A_clinical[which(KIF1A_clinical$Regression_hand == "Yes" | 
                     KIF1A_clinical$Regression_language == "Yes" | 
                     KIF1A_clinical$Regression_daily_living == "Yes" | 
                     KIF1A_clinical$Regression_motor == "Yes"), "Developmental_concern"] <- "Yes";
KIF1A_clinical[which(is.na(KIF1A_clinical$Developmental_concern)), "Developmental_concern"] <- "Unknown";

# Convert the date column to Date type
KIF1A_clinical$Date_survey <- as.Date(KIF1A_clinical$Date_survey);

# ------------------------- #
# KIF1A molecular phenotype #
# ------------------------- #

filename <- paste(rawdata_path, "/molecular_data/KIF1A mutants parameters (2026-04-02) v2.xlsx", sep="");
KIF1A_molecular <- read.xlsx(filename, sheet = 1);
KIF1A_molecular <- KIF1A_molecular[, setdiff(colnames(KIF1A_molecular), "conservation.old")];

colnames(KIF1A_molecular) <- c(
    "Variant_name", "Case_count", "Movement", "Diffusion", "Velocity", "Run_length", "Dwell_time", 
    "AA_ref", "AA_mut", "AA_pos", "Hydro_diff", "pI_diff", "MW_diff", "Dev_Z_score", "Weighted_penalty", 
    "Frac_ID", "Norm_ref_mean", "Norm_sop_mean", "Conservation");
rownames(KIF1A_molecular) <- KIF1A_molecular$Variant_name;
KIF1A_molecular$P_notation <- paste("p.", KIF1A_molecular$Variant_name, sep="");

# Remove variants not been purified (by Lu Rao)
KIF1A_molecular <- KIF1A_molecular[KIF1A_molecular$Variant_name != "K161M", ];

# =================== #
# Nextflow annotation #
# =================== #

filename <- paste(rawdata_path, "/predicted_score/KIF1A_VEP.vcf", sep="");
KIF1A.VEP <- read.table(filename, header = TRUE, stringsAsFactors = FALSE, sep = "\t", comment.char = "", skip = 5);
colnames(KIF1A.VEP) <- c("Chr", "Pos_start", "ID", "Ref_base", "Alt_base", "QUAL", "FILTER", "INFO");

filename <- paste(rawdata_path, "/predicted_score/KIF1A_variant_anno.tsv.gz", sep="");
if (file.exists(filename) & force_run == FALSE) {
    KIF1A.dnv_anno <- fread(file=filename, data.table=FALSE, quote=FALSE, sep="\t");
} else {
    vcf_file <- paste(rawdata_path, "/predicted_score/KIF1A.vcf", sep="");
    anno_path <- paste(rawdata_path, "/KIF1A_denovo_anno", sep="");
    KIF1A.dnv_anno <- de_novo_annotation(KIF1A.VEP, vcf_file, filename, nextflow_pars, anno_path);
}

var_names.raw <- setdiff(colnames(KIF1A.dnv_anno), c("ID", "CHROM", "POS", "REF", "ALT"));
var_names.raw <- KIF1A_anno[KIF1A_anno$Var_group == "Annotation" & KIF1A_anno$Var_raw_name %in% var_names.raw, ]$Var_raw_name;
var_names.new <- KIF1A_anno[KIF1A_anno$Var_group == "Annotation" & KIF1A_anno$Var_raw_name %in% var_names.raw, ]$Var_name;
KIF1A.dnv_anno <- KIF1A.dnv_anno[, c("ID", var_names.raw)];
colnames(KIF1A.dnv_anno) <- c("ID", var_names.new);

KIF1A.dnv_anno$P_notation <- KIF1A.dnv_anno$Canon_protein_change;
KIF1A.dnv_anno[grepl("[0-9]+[A-Z]>[0-9]+[A-Z]", KIF1A.dnv_anno$P_notation) == FALSE, ]$P_notation <- ".";
KIF1A.dnv_anno$P_notation <- gsub(">", "", gsub("^[0-9]+", "p.", KIF1A.dnv_anno$P_notation));

# Annotated MisFit v1.5 score (gnomAD version)
filename <- paste(rawdata_path, "/predicted_score/KIF1A_MisFit_gnomAD4.tsv", sep="");
KIF1A_MisFit <- read.table(file=filename, header=TRUE, stringsAsFactors=FALSE, quote=NULL, sep="\t");
KIF1A_MisFit$P_notation <- paste("p.", KIF1A_MisFit$AA_ref, KIF1A_MisFit$Protein_position, KIF1A_MisFit$AA_alt, sep="");
KIF1A_MisFit <- unique(KIF1A_MisFit[, c("P_notation", "MisFit_D", "MisFit_S")]);

# Replace MisFit v1.0 to MisFit v1.5
KIF1A.dnv_anno <- KIF1A.dnv_anno[, setdiff(colnames(KIF1A.dnv_anno), c("MisFit_D", "MisFit_S", "MisFit_S_gene"))];
KIF1A.dnv_anno <- merge(KIF1A.dnv_anno, KIF1A_MisFit, by = "P_notation", all.x = TRUE, all.y = FALSE);
KIF1A.dnv_anno$MisFit_S_gene <- 0.454;

# Annotated ESM score (https://huggingface.co/spaces/ntranoslab/esm_variants)
filename <- paste(rawdata_path, "/predicted_score/KIF1A_ESM.csv", sep="");
KIF1A_ESM <- read.table(file=filename, header=TRUE, stringsAsFactors=FALSE, quote=NULL, sep=",");
colnames(KIF1A_ESM) <- c("Variant", "ESM", "Position");
KIF1A_ESM$P_notation <- paste("p.", KIF1A_ESM$Variant, sep="");

# Map ESM score to the variant
KIF1A.dnv_anno <- merge(KIF1A.dnv_anno, KIF1A_ESM[, c("P_notation", "ESM")], by = "P_notation", all.x = TRUE, all.y = FALSE);

# ============================================================ #
# Merge KIF1A clinical data, variants data, and molecular data #
# ============================================================ #

# === Merge variant and clinical data ===
KIF1A_merge <- merge(KIF1A_variant, KIF1A_clinical, by = "Record_ID", all.x = TRUE, all.y = TRUE);

# === Merge with annotation data ===
KIF1A_merge <- merge(KIF1A_merge, KIF1A.dnv_anno, by = "ID", all.x = TRUE, all.y = FALSE);

# === Merge molecular data ===
KIF1A_merge <- merge(KIF1A_merge, KIF1A_molecular, by = "P_notation", all.x = TRUE, all.y = TRUE);
KIF1A_merge <- KIF1A_merge[!is.na(KIF1A_merge$Record_ID) & !is.na(KIF1A_merge$Instance), ];
KIF1A_merge$AA_pos <- as.numeric(gsub("[a-zA-Z.]+", "", KIF1A_merge$P_notation));

# Correct variant name
KIF1A_merge$Variant_name <- gsub("p.", "", KIF1A_merge$P_notation);
KIF1A_merge[which(KIF1A_merge$Variant_name == "."), ]$Variant_name <- NA;

# =================================== #
# Define new variables in KIF1A_merge #
# =================================== #

unique_IDs <- unique(KIF1A_merge[KIF1A_merge$Instance != 0, ]$Record_ID);
select_colnames <- KIF1A_anno[KIF1A_anno$Var_group == "Clinical" & KIF1A_anno$Var_type == 2, ]$Var_name;
select_colnames <- intersect(select_colnames, colnames(KIF1A_merge));

# ------------------------- #
# Define lifetime variables #
# ------------------------- #
new_lines <- list();
for (ID in unique_IDs) {
    temp_data <- KIF1A_merge[KIF1A_merge$Record_ID == ID, ];
    new_lines[[ID]] <- temp_data[which.max(temp_data$Instance), ];
    new_lines[[ID]]$Instance <- 0;
    for (name in select_colnames) {
        new_lines[[ID]][, name] <- "Unknown";
        temp_value <- temp_data[, name];
        temp_value <- temp_value[!is.na(temp_value)];
        if (length(temp_value) > 0) {
            if (any(temp_value == "Yes")) {
                new_lines[[ID]][, name] <- "Yes";
            } else if (any(temp_value == "Abnormal")) {
                new_lines[[ID]][, name] <- "Abnormal";
            } else if (any(temp_value == "Permanent")) {
                new_lines[[ID]][, name] <- "Permanent";
            } else if (all(temp_value == "No")) {
                new_lines[[ID]][, name] <- "No";
            } else if (all(temp_value == "Normal")) {
                new_lines[[ID]][, name] <- "Normal";
            } else if (all(temp_value == "Temporary")) {
                new_lines[[ID]][, name] <- "Temporary";
            }
        }
    }
}
KIF1A_merge <- rbind(KIF1A_merge, do.call(rbind, new_lines));
rm(new_lines);

# ------------------------------------------------------------- #
# Define final VABS, change of VABS, and rate of change of VABS #
# ------------------------------------------------------------- #

unique_IDs <- unique(KIF1A_merge[!is.na(KIF1A_merge$Age_VABS), ]$Record_ID);
VABS_vars <- KIF1A_anno[KIF1A_anno$Var_group == "VABS" & KIF1A_anno$Var_subgroup == "VABS", ]$Var_name;
VABS_diff_vars <- paste(VABS_vars, "_Diff", sep = "");
VABS_rate_vars <- paste(setdiff(VABS_vars, "Age_VABS"), "_Rate", sep = "");

# === Final VABS ===
for (ID in unique_IDs) {
    KIF1A_temp <- KIF1A_merge[KIF1A_merge$Record_ID %in% ID, ];
    KIF1A_temp <- KIF1A_temp[!is.na(KIF1A_temp$Age_VABS), ];
    KIF1A_temp <- KIF1A_temp[order(KIF1A_temp$Instance, decreasing = TRUE), ];
    temp_value <- as.numeric(KIF1A_temp[1, VABS_vars]);
    KIF1A_merge[KIF1A_merge$Record_ID == ID & KIF1A_merge$Instance == 0, VABS_vars] <- temp_value;
}

# === Change of VABS ===
for (name in VABS_diff_vars) {
    KIF1A_merge[[name]] <- NA;
}
for (ID in unique_IDs) {
    KIF1A_temp <- KIF1A_merge[!is.na(KIF1A_merge$Age_VABS) & 
                              KIF1A_merge$Record_ID == ID & 
                              KIF1A_merge$Instance > 0, ];
    if (nrow(KIF1A_temp) >= 2) {
        KIF1A_temp <- KIF1A_temp[order(KIF1A_temp$Instance, decreasing = FALSE), ];
        for (i in c(2:nrow(KIF1A_temp))) {
            temp_value <- as.numeric(KIF1A_temp[i, VABS_vars]) - as.numeric(KIF1A_temp[1, VABS_vars]);
            KIF1A_merge[KIF1A_merge$Record_ID == ID & KIF1A_merge$Instance == i, VABS_diff_vars] <- temp_value;
        }
        # Lifetime change of VABS
        temp_value <- KIF1A_merge[KIF1A_merge$Record_ID == ID & KIF1A_merge$Instance == nrow(KIF1A_temp), VABS_diff_vars];
        KIF1A_merge[KIF1A_merge$Record_ID == ID & KIF1A_merge$Instance == 0, VABS_diff_vars] <- as.numeric(temp_value);
    }
}

# === Rate of change of VABS ===
for (name in VABS_rate_vars) {
    KIF1A_merge[[name]] <- NA;
}
for (ID in unique_IDs) {
    KIF1A_temp <- KIF1A_merge[!is.na(KIF1A_merge$Age_VABS) & 
                              KIF1A_merge$Record_ID == ID & 
                              KIF1A_merge$Instance > 0, ];
    if (nrow(KIF1A_temp) >= 2) {
        KIF1A_temp <- KIF1A_temp[order(KIF1A_temp$Instance, decreasing = FALSE), ];
        for (i in c(2:nrow(KIF1A_temp))) {
            temp_value <- as.numeric(KIF1A_temp[i, setdiff(VABS_diff_vars, "Age_VABS_Diff")]) / 
                          as.numeric(KIF1A_temp[i, "Age_VABS_Diff"]);
            KIF1A_merge[KIF1A_merge$Record_ID == ID & KIF1A_merge$Instance == i, VABS_rate_vars] <- temp_value;
        }
        # Lifetime rate of change of VABS
        temp_value <- KIF1A_merge[KIF1A_merge$Record_ID == ID & KIF1A_merge$Instance == nrow(KIF1A_temp), VABS_rate_vars];
        KIF1A_merge[KIF1A_merge$Record_ID == ID & KIF1A_merge$Instance == 0, VABS_rate_vars] <- as.numeric(temp_value);
    }
}

# -------------------- #
# Define new variables #
# -------------------- #

KIF1A_merge$EEG_or_seizure <- "Unknown";
KIF1A_merge[KIF1A_merge$Seizure_any == "No" & KIF1A_merge$EEG_result == "Normal", ]$EEG_or_seizure <- "Normal";
KIF1A_merge[KIF1A_merge$Seizure_any == "Yes" | KIF1A_merge$EEG_result == "Abnormal", ]$EEG_or_seizure <- "Abnormal";

KIF1A_merge$Dwell_time_group <- KIF1A_merge$Dwell_time;
KIF1A_merge[!is.na(KIF1A_merge$Dwell_time) & KIF1A_merge$Dwell_time < 15, ]$Dwell_time_group <- "< 15";
KIF1A_merge[!is.na(KIF1A_merge$Dwell_time) & KIF1A_merge$Dwell_time >= 15, ]$Dwell_time_group <- "≥ 15";

# Only keep 1st MHI and last MHI for each cases
KIF1A_merge.endpoint <- do.call(rbind, lapply(split(
    KIF1A_merge,
    KIF1A_merge$Record_ID
), function(x) {
    # Find last Instance > 0 with non-NA VABS_ABC
    x_valid <- x[x$Instance > 0 & !is.na(x$VABS_ABC), ];
    if (nrow(x_valid) == 0) return(x[x$Instance == 0, ]);
    max_valid_instance <- max(x_valid$Instance);
    # Set intermediate instances (between first valid and last valid) to 0
    x$Instance[x$Instance > 1 & x$Instance < max_valid_instance] <- 0;
    # Also set instances after last valid to 0
    x$Instance[x$Instance > max_valid_instance] <- 0;
    return(x);
}));
KIF1A_merge.endpoint <- KIF1A_merge.endpoint[KIF1A_merge.endpoint$Instance > 0, ];
KIF1A_merge.endpoint[KIF1A_merge.endpoint$Instance > 1, ]$Instance <- 2;

# ================================== #
# Map the KIF1A_merge and KIF1A_anno #
# ================================== #

# === Rename columns ===
colnames(KIF1A_merge)[colnames(KIF1A_merge) == "ID"] <- "Variant_ID";
colnames(KIF1A_merge)[colnames(KIF1A_merge) == "K_ID"] <- "Historical_ID";
KIF1A_merge <- KIF1A_merge[, KIF1A_anno[KIF1A_anno$Var_name %in% colnames(KIF1A_merge), ]$Var_name];

KIF1A_merge.map <- map_pData_anno(KIF1A_merge, KIF1A_anno);
for (i in c(1:nrow(KIF1A_anno))) {
    if (KIF1A_anno$Var_type[i] == "1" | KIF1A_anno$Var_type[i] == "2") {
        var_name <- KIF1A_anno$Var_name[i];
        if (var_name %in% colnames(KIF1A_merge.map)) {
            KIF1A_merge.map[, var_name] <- as.numeric(KIF1A_merge.map[, var_name]);
        }
    }
}

# Calculate NA_count and NA_percentage in KIF1A_anno
for (i in c(1:nrow(KIF1A_anno))) {
    var_name <- KIF1A_anno$Var_name[i];
    if (var_name %in% colnames(KIF1A_merge.map)) {
        temp_val <- KIF1A_merge.map[, var_name];
        KIF1A_anno$NA_count[i] <- length(temp_val[is.na(temp_val)]);
    } else {
        KIF1A_anno$NA_count[i] <- nrow(KIF1A_merge.map);
    }
    KIF1A_anno$NA_percentage[i] <- KIF1A_anno$NA_count[i] / nrow(KIF1A_merge.map);
}

# Factor variables
var_list <- KIF1A_anno[KIF1A_anno$Var_group == "Clinical" & KIF1A_anno$Var_type == 2, ]$Var_name;
var_list <- intersect(var_list, colnames(KIF1A_merge));
for (var in var_list) {
    KIF1A_merge[, var] <- factor(KIF1A_merge[, var], levels = get_var_annotation(KIF1A_anno, var)$Detail);
}

# ==================== #
# Predefined Variables #
# ==================== #

clin_vars <- c(
    "Optic_nerve_change", "Seizure_any", "EEG_result", "Hypertonia", "Hypotonia", 
    "Microcephaly", "Ataxia", "Dystonia", "Neuropathy", "Developmental_concern"
);
patho_vars <- c(
    "REVEL", "gMVP", "ESM", "AlphaMissense", "MisFit_S", "MisFit_D"
);
mol_vars <- c(
    "Movement", "Diffusion", "Velocity", "Run_length", "Dwell_time", 
    "Frac_ID", "Norm_ref_mean", "Norm_sop_mean", "Conservation", 
    "Hydro_diff", "pI_diff", "MW_diff", "Dev_Z_score"
);

color_list <- list(
    tab10 = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", 
              "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF"), 
    tab20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C", 
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5", 
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F", 
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")
);

color.seq1 <- list(
    A = c("#FEE5D9", "#FCAE91", "#FB6A4A", "#DE2D26", "#A50F15"), 
    B = c("#FEEDDE", "#FDCC8A", "#FDAE6B", "#F7934E", "#D95F0E"), 
    C = c("#EFF3FF", "#BDD7E7", "#6BAED6", "#3182BD", "#08519C"), 
    D = c("#EDF8E9", "#BAE4B3", "#74C476", "#31A354", "#006D2C"), 
    E = c("#F2F0F7", "#CBC9E2", "#9E9AC8", "#756BB1", "#54278F")
);

color.seq2 <- list(
    A = c("#FEF0D9", "#FDCC8A", "#FC8D59", "#E34A33", "#B30000"), 
    B = c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404"), 
    C = c("#FFFFB2", "#FECC5C", "#FD8D3C", "#F03B20", "#BD0026"), 
    D = c("#FFFFCC", "#C2E699", "#78C679", "#31A354", "#006837"), 
    E = c("#EDF8FB", "#B2E2E2", "#66C2A4", "#2CA25F", "#006D2C"), 
    F = c("#F6EFF7", "#BDC9E1", "#67A9CF", "#1C9099", "#016C59"), 
    G = c("#F0F9E8", "#BAE4BC", "#7BCCC4", "#43A2CA", "#0868AC"), 
    H = c("#FFFFCC", "#A1DAB4", "#41B6C4", "#2C7FB8", "#253494"), 
    I = c("#F1EEF6", "#BDC9E1", "#74A9CF", "#2B8CBE", "#045A8D"), 
    J = c("#EDF8FB", "#B3CDE3", "#8C96C6", "#8856A7", "#810F7C"), 
    K = c("#F1EEF6", "#D7B5D8", "#DF65B0", "#DD1C77", "#980043"), 
    L = c("#FEEBE2", "#FBB4B9", "#F768A1", "#C51B8A", "#7A0177")
);

# Cases with all VABS ABC >= 100
removed_IDs <- c();
for (ID in unique(KIF1A_merge$Record_ID)) {
    temp_value <- KIF1A_merge[KIF1A_merge$Record_ID == ID, ]$VABS_ABC;
    temp_value <- temp_value[!is.na(temp_value)];
    if (length(temp_value) > 0) {
        if (all(temp_value >= 100)) {
            removed_IDs <- c(removed_IDs, ID);
        }
    }
}

# ========================== #
# Trajectory of VABS and age #
# ========================== #

VABS_trend_path <- paste(project_path, "/2_VABS_trend", sep="");
if (!dir.exists(VABS_trend_path)) {
    dir.create(VABS_trend_path);
}

KIF1A_temp <- KIF1A_merge.map[KIF1A_merge.map$Instance > 0, ];

# Overall trend of VABS_ABC by different variables
filepath <- paste(VABS_trend_path, "/Overall_trend", sep="");
if (!dir.exists(filepath)) {
    dir.create(filepath);
}

var_list <- c("Sex", "Death", "P_notation", "Movement", "Diffusion");
for (var_name in var_list) {
    message(var_name);
    filename <- paste(VABS_trend_path, "/VABS_ABC.", var_name, ".png", sep="");
    VABS_trend.line_plot(KIF1A_temp, KIF1A_anno, "VABS_ABC", var_name, color_list$tab10, filename = filename);
}

var_list <- c("REVEL", "ESM", "AlphaMissense", "MisFit_S", "MisFit_D");
for (i in c(1:length(var_list))) {
    var_name <- var_list[i];
    message(var_name);
    filename <- paste(VABS_trend_path, "/VABS_ABC.", var_name, ".png", sep="");
    VABS_trend.line_plot(KIF1A_temp, KIF1A_anno, "VABS_ABC", var_name, color.seq1[[i]], filename = filename);
}

var_list <- c("Velocity", "Run_length", "Dwell_time", "Conservation", "Hydro_diff", "pI_diff", "MW_diff");
for (i in c(1:length(var_list))) {
    var_name <- var_list[i];
    message(var_name);
    filename <- paste(VABS_trend_path, "/VABS_ABC.", var_name, ".png", sep="");
    VABS_trend.line_plot(KIF1A_temp, KIF1A_anno, "VABS_ABC", var_name, color.seq2[[i]], filename = filename);
}

# Separate by AA_pos
filepath <- paste(VABS_trend_path, "/Separate_by_AA_pos", sep="");
if (!dir.exists(filepath)) {
    dir.create(filepath);
}

unique_AA_pos <- unique(KIF1A_temp[!is.na(KIF1A_temp$VABS_ABC) & !is.na(KIF1A_temp$AA_pos), ]$AA_pos);
unique_AA_pos <- sort(unique_AA_pos[!is.na(unique_AA_pos)]);
for (AA_pos in unique_AA_pos) {
    message(paste("AA_pos ", AA_pos, sep=""));
    filename <- paste(filepath, "/VABS_ABC.AA_pos_", AA_pos, ".png", sep="");
    temp_data <- KIF1A_temp[!is.na(KIF1A_temp$VABS_ABC) & !is.na(KIF1A_temp$AA_pos), ];
    VABS_trend.line_plot(temp_data[temp_data$AA_pos == AA_pos, ], 
        KIF1A_anno, "VABS_ABC", "P_notation", color_list$tab10, c(0, 60), c(20, 140), TRUE, 3.8, 6, filename);
}

# ------------------------- #
# Highly recurrent variants #
# ------------------------- #

# === Variant frequency ===
variant_freq <- as.data.frame(table(KIF1A_merge[KIF1A_merge$Instance == 0, "P_notation"], dnn = "P_notation"));
variant_freq <- variant_freq[variant_freq$P_notation != ".", ];

temp_data <- KIF1A_merge[KIF1A_merge$Instance == 0, c("P_notation", "Variant_name", "REVEL", "ESM", "MisFit_D", "MisFit_S")];
temp_data <- unique(temp_data[!is.na(temp_data$P_notation), ]);

variant_freq <- merge(variant_freq, temp_data, by = "P_notation", all.x = TRUE, all.y = FALSE);
variant_freq <- variant_freq[order(variant_freq$Freq, decreasing = TRUE), ];

# Highly recurrent variants
filepath <- paste(VABS_trend_path, "/Highly_recurrent_variants", sep="");
if (!dir.exists(filepath)) {
    dir.create(filepath);
}

filename <- paste(filepath, "/VABS_ABC.min_freq_5.png", sep="");
VABS_trend.scatter_plot(KIF1A_merge, variant_freq, min_freq = 5, colors_list = color_list$tab20, 
                        show_legend = TRUE, height = 5, width = 7.5, filename = filename);

filename <- paste(filepath, "/VABS_ABC.min_freq_5.no_legend.png", sep="");
VABS_trend.scatter_plot(KIF1A_merge, variant_freq, min_freq = 5, colors_list = color_list$tab20, 
                        show_legend = FALSE, height = 5, width = 7, filename = filename);

# ---------------------- #
# Clustering by patients #
# ---------------------- #

filepath <- paste(VABS_trend_path, "/Cluster_by_patient.k_2", sep="");
patient_cluster <- VABS_trend.cluster(KIF1A_merge, vabs_col = "VABS_ABC", k = 2, mean_vabs_weight = 3, slope_weight = 1, 
                                      color_list = color_list$tab10[c(4, 1, 3)], filepath = filepath);

filepath <- paste(VABS_trend_path, "/Cluster_by_patient.k_3", sep="");
patient_cluster <- VABS_trend.cluster(KIF1A_merge, vabs_col = "VABS_ABC", k = 3, mean_vabs_weight = 3, slope_weight = 1, 
                                      color_list = color_list$tab10[c(4, 1, 3)], filepath = filepath);

filepath <- paste(VABS_trend_path, "/Cluster_by_patient.endpoint.k_2", sep="");
patient_cluster <- VABS_trend.cluster(KIF1A_merge.endpoint, vabs_col = "VABS_ABC", k = 2, mean_vabs_weight = 3, slope_weight = 1, 
                                      color_list = color_list$tab10[c(4, 1, 3)], filepath = filepath);

filepath <- paste(VABS_trend_path, "/Cluster_by_patient.endpoint.k_3", sep="");
patient_cluster <- VABS_trend.cluster(KIF1A_merge.endpoint, vabs_col = "VABS_ABC", k = 3, mean_vabs_weight = 3, slope_weight = 1, 
                                      color_list = color_list$tab10[c(4, 1, 3)], filepath = filepath);

# Add Patient_cluster to KIF1A_merge
temp_data <- data.frame(
    Record_ID       = patient_cluster$patient_features$Record_ID,
    Patient_cluster = patient_cluster$patient_features$cluster, 
    stringsAsFactors = FALSE
);
if (!("Patient_cluster" %in% colnames(KIF1A_merge))) {
    KIF1A_merge <- merge(KIF1A_merge, temp_data, by = "Record_ID", all.x = TRUE, all.y = FALSE);
}

# Separate by seizure resistant or not
temp_data <- KIF1A_merge.map[KIF1A_merge.map$Instance == 0 & KIF1A_merge.map$Seizure_resistant != 1, ];
temp_data$Seizure_resistant_group <- temp_data$Seizure_resistant;
temp_data[temp_data$Seizure_resistant_group == 2, "Seizure_resistant_group"] <- 1;

KIF1A_temp <- KIF1A_merge.map[KIF1A_merge.map$Instance != 0 & KIF1A_merge.map$Record_ID %in% temp_data$Record_ID, ];
KIF1A_temp <- merge(KIF1A_temp, temp_data[, c("Record_ID", "Seizure_resistant_group")], 
                    by = "Record_ID", all.x = TRUE, all.y = FALSE);

filename <- paste(VABS_trend_path, "/VABS_ABC_by_Seizure_resistant.png", sep="");
VABS_trend.line_plot(KIF1A_temp, KIF1A_anno, "VABS_ABC", "Seizure_resistant_group", color_list$tab10, filename = filename);

# ------------------------------- #
# Group-based trajectory modeling #
# ------------------------------- #

filepath <- file.path(VABS_trend_path, "GBTM_cluster");
VABS_trend.all <- VABS_trend.GBTM(
    KIF1A_merge,
    vabs_col   = "VABS_ABC",
    ng_max     = 4,
    log_age    = TRUE,
    color_list = color_list$tab10[c(2, 10, 3, 1)],
    filepath   = filepath
);

temp_IDs <- KIF1A_merge[KIF1A_merge$Instance == 2 & !is.na(KIF1A_merge$Age_VABS), ]$Record_ID;
KIF1A_temp <- KIF1A_merge[(KIF1A_merge$Record_ID %in% temp_IDs & KIF1A_merge$Instance < 2) | KIF1A_merge$Instance == 1, ];
KIF1A_temp[KIF1A_temp$Instance == 0, "Instance"] <- 2;

filepath <- file.path(VABS_trend_path, "GBTM_cluster.lifetime");
VABS_trend.lifetime <- VABS_trend.GBTM(
    KIF1A_temp,
    vabs_col   = "VABS_ABC",
    ng_max     = 4,
    log_age    = TRUE,
    color_list = color_list$tab10[c(2, 10, 3, 1)],
    filepath   = filepath
);

# === Add GBTM cluster results to KIF1A_merge and KIF1A_merge.map ===
temp_data <- VABS_trend.all[[2]]$pprob[, c("Record_ID", "class")];
colnames(temp_data) <- c("Record_ID", "GBTM_cluster");
if (!("GBTM_cluster" %in% colnames(KIF1A_merge))) {
    KIF1A_merge <- merge(KIF1A_merge, temp_data, by = "Record_ID", all.x = TRUE, all.y = FALSE);
}
if (!("GBTM_cluster" %in% colnames(KIF1A_merge.map))) {
    KIF1A_merge.map <- merge(KIF1A_merge.map, temp_data, by = "Record_ID", all.x = TRUE, all.y = FALSE);
}

# =================== #
# Clustering analysis #
# =================== #

cluster_path <- file.path(project_path, "3_Clustering_analysis");
if (!dir.exists(cluster_path)) {
    dir.create(cluster_path);
}

KIF1A_molecular.scaled <- KIF1A_molecular[, mol_vars];
KIF1A_molecular.scaled <- KIF1A_molecular.scaled[which(rownames(KIF1A_molecular.scaled) != "WT"), ];

KIF1A_molecular.scaled <- apply(KIF1A_molecular.scaled, 2, function(x) {
    (x - min(x)) / (max(x) - min(x))
});
KIF1A_molecular.scaled[, c(3:13)] <- signif(KIF1A_molecular.scaled[, c(3:13)], 6);

filename <- file.path(cluster_path, "KIF1A_molecular_scaled_data.txt");
write.table(KIF1A_molecular.scaled, filename, quote=FALSE, sep="\t");

# ------------------------------------ #
# Heatmap of the scaled molecular data #
# ------------------------------------ #

filepath <- file.path(cluster_path, "Hierarchical_clustering");
if (!dir.exists(filepath)) {
    dir.create(filepath);
}

KIF1A_temp <- KIF1A_merge[KIF1A_merge$Instance == 0 & KIF1A_merge$P_notation %in% KIF1A_molecular$P_notation, ];

filename <- file.path(filepath, "KIF1A_molecular_scaled.heatmap.png");
molecular_ht.h <- molecular_data.heatmap(KIF1A_molecular.scaled, KIF1A_temp, "h", 6, 16, filename);

filename <- file.path(filepath, "KIF1A_molecular_scaled.heatmap.new.png");
molecular_ht.v <- molecular_data.heatmap(KIF1A_molecular.scaled, KIF1A_temp, "v", 12, 6, filename);

# === Add cluster results to KIF1A_merge and KIF1A_merge.map ===
clust_groups <- cutree(molecular_ht.h$hclust, k = 4);
temp_data <- data.frame(
    Variant_name      = names(clust_groups),
    Molecular_cluster = clust_groups, 
    stringsAsFactors = FALSE
);

# Remap cluster labels: 1->3, 2->1, 3->4, 4->2
cluster_map <- c("1" = 3, "2" = 1, "3" = 4, "4" = 2);
temp_data$Molecular_cluster <- cluster_map[as.character(temp_data$Molecular_cluster)];

if (!("Molecular_cluster" %in% colnames(KIF1A_merge))) {
    KIF1A_merge <- merge(KIF1A_merge, temp_data, by = "Variant_name", all.x = TRUE, all.y = FALSE);
}
if (!("Molecular_cluster" %in% colnames(KIF1A_merge.map))) {
    KIF1A_merge.map <- merge(KIF1A_merge.map, temp_data, by = "Variant_name", all.x = TRUE, all.y = FALSE);
}

# === Correlation between molecular features and pathogenicity scores ===
KIF1A_temp <- KIF1A_merge[KIF1A_merge$Instance == 0 & KIF1A_merge$P_notation %in% KIF1A_molecular$P_notation, ];
temp_data <- unique(KIF1A_temp[, c("Variant_name", "Molecular_cluster", patho_vars, mol_vars)]);

pairwise_heatmap <- function(
        matrix,
        var_list_1    = NULL,
        var_list_2    = NULL,
        color_list    = NULL,
        color_range   = NULL,
        sig_color     = "#FF0000",
        font_size     = 10,
        decimal       = 2,
        col_names_rot = 30,
        width         = 8,
        height        = 6,
        filename      = NULL
    ) {

    if (is.null(var_list_1)) var_list_1 <- colnames(matrix);
    if (is.null(var_list_2)) var_list_2 <- colnames(matrix);
    if (is.null(color_list)) color_list <- rev(brewer.pal(n = 7, name = "RdYlBu"));
    if (is.null(color_range)) color_range <- c(-1, 1);

    # ── Compute pairwise correlation and p-value matrices ─────────────────────
    cor_mat <- matrix(NA, nrow = length(var_list_1), ncol = length(var_list_2),
                      dimnames = list(var_list_1, var_list_2));
    p_mat   <- cor_mat;

    for (v1 in var_list_1) {
        for (v2 in var_list_2) {
            x <- as.numeric(matrix[, v1]);
            y <- as.numeric(matrix[, v2]);
            valid <- !is.na(x) & !is.na(y);
            if (sum(valid) >= 3) {
                res             <- cor.test(x[valid], y[valid], method = "pearson");
                cor_mat[v1, v2] <- res$estimate;
                p_mat[v1, v2]   <- res$p.value;
            };
        };
    };

    # ── Fill NA values for display ────────────────────────────────────────────
    na_mat          <- is.na(cor_mat);
    cor_mat[na_mat] <- 0;
    p_mat[na_mat]   <- 1;

    # ── Build cell labels and significance matrix ─────────────────────────────
    fmt <- paste0("%.", decimal, "f");

    label_mat <- matrix(
        sprintf(fmt, cor_mat),
        nrow = length(var_list_1),
        dimnames = list(var_list_1, var_list_2)
    );
    label_mat[na_mat] <- NA;

    sig_mat <- p_mat < 0.05;
    sig_mat[is.na(sig_mat)] <- FALSE;

    # ── ComplexHeatmap cell function ──────────────────────────────────────────
    cell_fun <- function(j, i, x, y, width, height, fill) {
        label <- label_mat[i, j];
        if (!is.na(label)) {
            if (sig_mat[i, j]) {
                grid::grid.text(label, x, y,
                                gp = grid::gpar(fontsize = font_size, col = sig_color));
            } else {
                grid::grid.text(label, x, y,
                                gp = grid::gpar(fontsize = font_size, col = "#A0A0A0"));
            };
        };
    };

    # ── Color bar breaks and labels ───────────────────────────────────────────
    color_breaks  <- seq(color_range[1], color_range[2], length.out = 100);
    legend_at     <- seq(color_range[1], color_range[2], length.out = 5);
    legend_labels <- sprintf(fmt, legend_at);

    # ── Draw heatmap ──────────────────────────────────────────────────────────
    ht <- ComplexHeatmap::Heatmap(
        cor_mat,
        col              = circlize::colorRamp2(color_breaks, colorRampPalette(color_list)(100)),
        cell_fun         = cell_fun,
        cluster_rows     = FALSE,
        cluster_columns  = FALSE,
        row_names_side   = "left",
        name             = "Pearson r",
        rect_gp          = grid::gpar(col = "white", lwd = 0.5),
        row_names_gp     = grid::gpar(fontsize = font_size),
        column_names_gp  = grid::gpar(fontsize = font_size),
        column_names_rot = col_names_rot,
        heatmap_legend_param = list(
            title         = "Pearson r",
            at            = legend_at,
            labels        = legend_labels,
            legend_height = grid::unit(4, "cm")
        )
    );

    # ── Output ────────────────────────────────────────────────────────────────
    if (is.null(filename)) {
        ComplexHeatmap::draw(ht);
    } else {
        ext <- tolower(tools::file_ext(filename));
        if (ext == "pdf") {
            pdf(filename, width = width, height = height);
        } else {
            png(filename, width = width, height = height, units = "in", res = 300);
        };
        ComplexHeatmap::draw(ht);
        dev.off();
    };

    invisible(list(cor = cor_mat, p = p_mat, na = na_mat));
};

filename <- file.path(filepath, "Heatmap.patho_mol.overall.png");
pairwise_heatmap(temp_data, c(patho_vars, mol_vars), c(patho_vars, mol_vars), 
                 color_range = c(-1.0, 1.0), height = 8.5, width = 10, filename = filename);

filename <- file.path(filepath, "Heatmap.patho_mol.pairwise.png");
pairwise_heatmap(temp_data, mol_vars, patho_vars, color_list = c("#809EC2", "#FFFFFF", "#F3A447"), 
                 color_range = c(-0.6, 0.6), height = 5, width = 5, filename = filename);

for (i in unique(temp_data$Molecular_cluster)) {
    filename <- file.path(filepath, paste("Heatmap.patho_mol.pairwise.cluster_", i, ".png", sep=""));
    pairwise_heatmap(temp_data[temp_data$Molecular_cluster == i, ], mol_vars, patho_vars, 
                     color_list = c("#809EC2", "#FFFFFF", "#F3A447"), 
                     color_range =c(-0.6, 0.6), height = 5, width = 5, filename = filename);
}

# ---------------------------------- #
# PCA of molecular and clinical data #
# ---------------------------------- #

# === Molecular data ===
out_path <- file.path(cluster_path, "PCA.molecular.k_2");
PCA_reslut_1 <- PCA_analysis(KIF1A_molecular.scaled, out_path, 8, 2);

out_path <- file.path(cluster_path, "PCA.molecular.k_4");
PCA_reslut_1 <- PCA_analysis(KIF1A_molecular.scaled, out_path, 8, 4);

# === Molecular + clinical data ===
select_colnames <- c("Sex", "Age_VABS", "VABS_ABC", clin_vars, mol_vars);
temp_data <- KIF1A_merge.map[KIF1A_merge.map$Instance == 0 & !is.na(KIF1A_merge.map$Variant_name), ];
temp_data <- temp_data[, c("Record_ID", select_colnames)];
rownames(temp_data) <- temp_data$Record_ID;
temp_data <- temp_data[, setdiff(colnames(temp_data), "Record_ID")];
temp_data <- temp_data[!is.na(apply(temp_data, 1, mean)), ];

out_path <- file.path(cluster_path, "PCA.molecular_clincial.k_2");
PCA_reslut_2 <- PCA_analysis(temp_data, out_path, 8, 2);

out_path <- file.path(cluster_path, "PCA.molecular_clincial.k_4");
PCA_reslut_2 <- PCA_analysis(temp_data, out_path, 8, 4);

# === Save PCA clusters to file ===
temp_data <- KIF1A_merge[KIF1A_merge$Instance == 0, c("Record_ID", "Variant_name", select_colnames)];
temp_data <- merge(temp_data, PCA_reslut_1$result[, c("ID", "PCA_cluster")], 
                   by.x = "Variant_name", by.y = "ID", all.x = TRUE, all.y = FALSE);
colnames(temp_data)[colnames(temp_data) == "PCA_cluster"] <- "PCA_cluster.molecular";
temp_data <- merge(temp_data, PCA_reslut_2$result[, c("ID", "PCA_cluster")], 
                   by.x = "Record_ID", by.y = "ID", all.x = TRUE, all.y = FALSE);
colnames(temp_data)[colnames(temp_data) == "PCA_cluster"] <- "PCA_cluster.molecular_clinical";

write.xlsx(temp_data, file.path(cluster_path, "KIF1A_PCA_cluster.xlsx"), rowNames = FALSE);

# ============ #
# Scatter plot #
# ============ #

scatter_plot_path <- paste(project_path, "/4_Scatter_plot", sep="");
if (!dir.exists(scatter_plot_path)) {
    dir.create(scatter_plot_path);
}

scatter_plot <- function(data, x, y, group = NULL, x_lim = NULL, y_lim = NULL,
                         x_lab = NULL, y_lab = NULL, x_log = FALSE, y_log = FALSE,
                         pt_size = 2.5, pt_alpha = 0.7, pt_color = "#4F81BD",
                         lm_fit = FALSE, lm_color = "#C0504D", show_P_val = TRUE,
                         height = 4.5, width = 6, filename = NULL) {
    # Remove NA rows
    keep_cols <- c(x, y, if (!is.null(group)) group);
    data <- data[complete.cases(data[, keep_cols]), ];

    # Build base plot
    if (!is.null(group)) {
        group_count <- length(unique(data[, group]));
        if (group_count > length(pt_color)) {
            pt_color <- rainbow(group_count);
        }
        p <- ggplot(data, aes(x = .data[[x]], y = .data[[y]], color = .data[[group]])) +
                theme_classic() +
                geom_point(size = pt_size, alpha = pt_alpha) +
                scale_color_manual(values = pt_color);
    } else {
        p <- ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
                theme_classic() +
                geom_point(size = pt_size, alpha = pt_alpha, color = pt_color);
    }

    # Log scale
    if (x_log) p <- p + scale_x_log10();
    if (y_log) p <- p + scale_y_log10();

    # Linear fit
    if (lm_fit) {
        p <- p + geom_smooth(formula = y ~ x, method = "lm", col = lm_color);
        # P value annotation
        if (show_P_val) {
            p <- p + stat_poly_eq(formula = y ~ x, color = "#FF0000",
                                aes(label = get_rr_p_label(after_stat(r.squared), after_stat(p.value))),
                                parse = TRUE, label.x = "left", label.y = "top");
        }
    }

    # Axis labels
    x_lab <- if (is.null(x_lab)) x else x_lab;
    y_lab <- if (is.null(y_lab)) y else y_lab;

    # Axis limits
    if (!is.null(x_lim)) p <- p + xlim(x_lim);
    if (!is.null(y_lim)) p <- p + ylim(y_lim);

    p <- p + labs(x = x_lab, y = y_lab) +
            theme(axis.title = element_text(color = "black"),
                  axis.text  = element_text(color = "black"),
                  axis.ticks = element_line(color = "black"));

    # Save or return
    if (!is.null(filename)) {
        ggsave(filename, p, height = height, width = width, units = "in");
    } else {
        return(p);
    }
}

# ---------------------------- #
# Scatter plot of VABS and age #
# ---------------------------- #

# Extract age at 1st evaluation
temp_data <- KIF1A_merge[KIF1A_merge$Instance == 1, c("Record_ID", "Age_VABS")];
colnames(temp_data) <- c("Record_ID", "Age_VABS_1st");

# Merge final test data and 1st MHI age
KIF1A_temp <- KIF1A_merge[KIF1A_merge$Instance == 0 & !is.na(KIF1A_merge$Age_VABS_Diff), ];
KIF1A_temp <- merge(KIF1A_temp, temp_data, by = "Record_ID", all.x = TRUE, all.y = FALSE);

filepath <- paste(scatter_plot_path, "/Age_vs_VABS", sep="");
if (!dir.exists(filepath)) {
    dir.create(filepath);
}

# Age at 1st evaluation vs. Change of VABS ABC
filename <- paste(filepath, "/Age_1st_MHI_vs_VABS_ABC_Diff.png", sep="");
p <- scatter_plot(KIF1A_temp, "Age_VABS_1st", "VABS_ABC_Diff", x_lab = "Age at 1st evaluation", 
                  y_lab = "Change of VABS ABC", x_log = TRUE, pt_size = 1.5) + 
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") + 
        geom_vline(xintercept = 10, linetype = "dashed", color = "red");
ggsave(filename, p, height = 3, width = 3, units = "in");

# Age at 1st evaluation vs. Rate of Change of VABS ABC
filename <- paste(filepath, "/Age_1st_MHI_vs_VABS_ABC_Rate.png", sep="");
p <- scatter_plot(KIF1A_temp, "Age_VABS_1st", "VABS_ABC_Rate", x_lab = "Age at 1st evaluation", 
                  y_lab = "Rate of change of VABS ABC", x_log = TRUE, pt_size = 1.5) + 
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") + 
        geom_vline(xintercept = 10, linetype = "dashed", color = "red");
ggsave(filename, p, height = 3, width = 3, units = "in");

# ------------------------------------------------------ #
# Scatter plot of molecular data and pathogenicity score #
# ------------------------------------------------------ #

select_colnames <- c("Variant_name", patho_vars, mol_vars, "Molecular_cluster");
KIF1A_temp <- KIF1A_merge[KIF1A_merge$Instance == 0, select_colnames];
KIF1A_temp <- unique(KIF1A_temp[!is.na(KIF1A_temp$Molecular_cluster), ]);

# === Overall correlation ===
filepath <- file.path(scatter_plot_path, "Pathogenicity_score_vs_molecular");
if (!dir.exists(filepath)) {
    dir.create(filepath);
}
for (var_1 in patho_vars) {
    for (var_2 in mol_vars) {
        p <- scatter_plot(KIF1A_temp, var_1, var_2, x_lab = gsub("_", "-", var_1), 
                          y_lab = get_var_description(KIF1A_anno, var_2), 
                          lm_fit = TRUE, pt_size = 1.5);
        filename <- file.path(filepath, paste(var_1, "_vs_", var_2, ".png", sep = ""));
        ggsave(filename, p, height = 3, width = 3, units = "in");
    }
}

# === Each cluster ===
for (i in unique(KIF1A_temp$Molecular_cluster)) {
    temp_data <- KIF1A_temp[KIF1A_temp$Molecular_cluster == i, ];
    filepath <- file.path(scatter_plot_path, paste("Pathogenicity_score_vs_molecular.cluster_", i, sep=""));
    if (!dir.exists(filepath)) {
        dir.create(filepath);
    }
    if (i <= 2) {
        var_list <- setdiff(mol_vars, c("Movement", "Diffusion"));
    } else if (i == 3) {
        var_list <- c("Dwell_time", "Hydro_diff", "pI_diff", "MW_diff");
    } else {
        var_list <- c("Dwell_time", "Conservation", "Hydro_diff", "pI_diff", "MW_diff");
    }
    for (var_1 in patho_vars) {
        for (var_2 in var_list) {
            p <- scatter_plot(temp_data, var_1, var_2, x_lab = gsub("_", "-", var_1), 
                            y_lab = get_var_description(KIF1A_anno, var_2), 
                            lm_fit = TRUE, pt_size = 1.5);
            filename <- file.path(filepath, paste(var_1, "_vs_", var_2, ".png", sep = ""));
            ggsave(filename, p, height = 3, width = 3, units = "in");
        }
    }
}


scatter_plot_temp <- function(data, var) {
    p <- ggplot(data[!is.na(data[, var]) & !is.na(data[, var_name]) & !is.infinite(data[, var]), ], 
           aes(x = get(var), y = get(var_name))) + theme_classic() + 
           geom_point(color = "#4F81BD", alpha = 0.5) + 
           geom_smooth(formula = y ~ x, method = "lm", col = "#C0504D") + 
           stat_poly_eq(formula = y ~ x, color = "#FF0000", 
                        aes(label = get_rr_p_label(..r.squared.., ..p.value.., length(var_list))), 
                        parse = TRUE, label.x = "left", label.y = "top");
    if (var_anno) {
        p <- p + xlab(get_var_description(KIF1A_anno, var)) + 
                 ylab(get_var_description(KIF1A_anno, var_name));
    } else {
        p <- p + xlab(var) + ylab(var_name);
    }
    return(p);
}

# ------------------------------------------------ #
# Scatter plot of VABS ABC and different variables #
# ------------------------------------------------ #

KIF1A_temp <- KIF1A_merge[KIF1A_merge$Instance == 0, ];

# Scatter plot of final VABS ABC and pathogenicity score
var_name <- "VABS_ABC";
var_anno <- FALSE;
var_list <- c("REVEL", "gMVP", "ESM", "AlphaMissense", "MisFit_S", "MisFit_D");

filename <- paste(scatter_plot_path, "/VABS_ABC_Final_vs_pathogenicity_score.png", sep="");
p_list <- lapply(var_list, scatter_plot_temp, data = KIF1A_temp);
p_code <- paste("ggarrange(", paste(paste("p_list[[", 1:length(p_list), "]]", sep=""), collapse=", "), ", ncol=3, nrow=2)", sep="");
p_combine <- eval(parse(text=p_code));
ggsave(filename, p_combine, width=9, height=5, units="in");

# Scatter plot of VABS_ABC_Final and molecular data
var_name <- "VABS_ABC";
var_anno <- FALSE;
var_list <- c("Movement", "Diffusion", "Velocity", "Run_length", "Dwell_time", 
              "Conservation", "Hydro_diff", "pI_diff", "MW_diff");

filename <- paste(scatter_plot_path, "/VABS_ABC_Final_vs_molecular_data.png", sep="");
p_list <- lapply(var_list, scatter_plot_temp, data = KIF1A_merge);
p_code <- paste("ggarrange(", paste(paste("p_list[[", 1:length(p_list), "]]", sep=""), collapse=", "), ", ncol=3, nrow=3)", sep="");
p_combine <- eval(parse(text=p_code));
ggsave(filename, p_combine, width=9, height=7.5, units="in");

# Scatter plot of VABS_ABC_Diff and molecular data
var_name <- "VABS_ABC_Diff";
var_list <- c("Movement", "Diffusion", "Velocity", "Run_length", "Dwell_time", 
              "Conservation", "Hydro_diff", "pI_diff", "MW_diff");

filename <- paste(scatter_plot_path, "/VABS_ABC_Diff_vs_molecular_data.png", sep="");
p_list <- lapply(var_list, scatter_plot_temp, data = KIF1A_merge);
p_code <- paste("ggarrange(", paste(paste("p_list[[", 1:length(p_list), "]]", sep=""), collapse=", "), ", ncol=3, nrow=3)", sep="");
p_combine <- eval(parse(text=p_code));
ggsave(filename, p_combine, width=9, height=7.5, units="in");

# === Scatter plot of VABS_ABC_Final and molecular data with Age_VABS_Final > 7 ===
var_name <- "VABS_ABC_Final";
var_list <- c("Movement", "Diffusion", "Velocity", "Run_length", "Dwell_time", 
              "Conservation", "Hydro_diff", "pI_diff", "MW_diff");

filename <- paste(scatter_plot_path, "/VABS_ABC_Final_vs_molecular_data.age_above_7.png", sep="");
p_list <- lapply(var_list, scatter_plot_temp, data = KIF1A_merge[!is.na(KIF1A_merge$Age_VABS_Final) & KIF1A_merge$Age_VABS_Final > 7, ]);
p_code <- paste("ggarrange(", paste(paste("p_list[[", 1:length(p_list), "]]", sep=""), collapse=", "), ", ncol=3, nrow=3)", sep="");
p_combine <- eval(parse(text=p_code));
ggsave(filename, p_combine, width=9, height=7.5, units="in");

# === Grouped by Movement & Diffusion ===
var_name <- "VABS_ABC_Final";
var_list <- c("REVEL", "ESM", "AlphaMissense", "MisFit_S", "MisFit_D", "Velocity", 
    "Run_length", "Dwell_time", "Conservation", "Hydro_diff", "pI_diff", "MW_diff");

# Grouped by Movement
for (i in c(0, 1)) {
    filename <- paste(scatter_plot_path, "/VABS_ABC_Final_vs_others.Movement_", i, ".png", sep="");
    p_list <- lapply(var_list, scatter_plot_temp, data = KIF1A_merge[!is.na(KIF1A_merge$Movement) & KIF1A_merge$Movement == i, ]);
    p_code <- paste("ggarrange(", paste(paste("p_list[[", 1:length(p_list), "]]", sep=""), collapse=", "), ", ncol=4, nrow=3)", sep="");
    p_combine <- eval(parse(text=p_code));
    ggsave(filename, p_combine, width=12, height=7.5, units="in");
}

# Grouped by Diffusion
for (i in c(0, 1)) {
    filename <- paste(scatter_plot_path, "/VABS_ABC_Final_vs_others.Diffusion_", i, ".png", sep="");
    p_list <- lapply(var_list, scatter_plot_temp, data = KIF1A_merge[!is.na(KIF1A_merge$Diffusion) & KIF1A_merge$Diffusion == i, ]);
    p_code <- paste("ggarrange(", paste(paste("p_list[[", 1:length(p_list), "]]", sep=""), collapse=", "), ", ncol=4, nrow=3)", sep="");
    p_combine <- eval(parse(text=p_code));
    ggsave(filename, p_combine, width=12, height=7.5, units="in");
}

# -------------------------------------------------------- #
# Scatter plot of 1st MHI VABS_ABC and Age_VABS by Cluster #
# -------------------------------------------------------- #
scatter_plot_temp <- function(data, cluster) {
    temp_data <- data[!is.na(data$Cluster) & data$Cluster == cluster, ];
    ggplot(temp_data, aes(x=Age_VABS, y=VABS_ABC)) + theme_classic() + 
           geom_point() + scale_x_log10() + 
           ggtitle(paste("Cluster ", cluster, sep="")) + 
           xlab("Age of 1st VABS eval") + ylab("VABS ABC, 1st MHI") + 
           geom_vline(xintercept = 7, linetype = "dashed", color = "red");
}

filename <- paste(scatter_plot_path, "/VABS_ABC_vs_Age_VABS.1st_MHI.Cluster.png", sep="");
p_list <- lapply(c(1, 2, 3, 4), scatter_plot_temp, data = KIF1A_merge[KIF1A_merge$Instance == 1, ]);
p_code <- paste("ggarrange(", paste(paste("p_list[[", 1:length(p_list), "]]", sep=""), collapse=", "), ", ncol=2, nrow=2)", sep="");
p_combine <- eval(parse(text=p_code));
ggsave(filename, p_combine, width=10, height=7.5, units="in");

# ------------------------------------------------------------------------- #
# Scatter plot of final VABS ABC and different variables grouped by Cluster #
# ------------------------------------------------------------------------- #

scatter_plot_temp <- function(data, cluster) {
    temp_data <- data[!is.na(data$Cluster) & data$Cluster == cluster, ];
    ggplot(temp_data, aes(x = get(var_name), y = VABS_ABC)) + theme_classic() + 
           geom_point(color = "#4F81BD", alpha = 0.5) + 
           ggtitle(paste("Cluster ", cluster, sep="")) + 
           xlab(get_var_description(KIF1A_anno, var_name)) + ylab("Final VABS ABC") + 
           geom_smooth(formula = y ~ x, method = "lm", col = "#C0504D") + 
           stat_poly_eq(formula = y ~ x, color = "#FF0000", 
                        aes(label = get_rr_p_label(..r.squared.., ..p.value.., length(var_list))), 
                        parse = TRUE, label.x = "left", label.y = "top");
}

var_list <- c("REVEL", "ESM", "MisFit_S", "MisFit_D", "Velocity", "Run_length", "Dwell_time");
for (var_name in var_list) {
    filename <- paste(scatter_plot_path, "/VABS_ABC_Final_vs_", var_name, ".Cluster.png", sep="");
    p_list <- lapply(c(1, 2, 3, 4), scatter_plot_temp, data = KIF1A_merge[KIF1A_merge$Instance == 0, ]);
    p_code <- paste("ggarrange(", paste(paste("p_list[[", 1:length(p_list), "]]", sep=""), collapse=", "), ", ncol=2, nrow=2)", sep="");
    p_combine <- eval(parse(text=p_code));
    ggsave(filename, p_combine, width=9, height=7.5, units="in");
}

# --------------------- #
# Specific scatter plot #
# --------------------- #

temp_data <- KIF1A_merge[, c("Run_length", "Dwell_time", "Variant_name")];
temp_data <- temp_data[!is.na(temp_data$Run_length) & !is.na(temp_data$Dwell_time), ];
temp_data <- temp_data[!duplicated(temp_data[, c("Run_length", "Dwell_time")]), ];
temp_data[temp_data$Run_length == 0 & temp_data$Dwell_time < 3, ]$Variant_name <- "";

color_values <- c("Abnormal" = "#E74C3C", "Normal/Unknown" = "#3498DB")

p1 <- ggplot(KIF1A_merge, aes(x = Run_length, y = Dwell_time)) + 
    theme_classic() + 
    geom_jitter(aes(color = EEG_seizure_lifetime), 
                width = 0.05, height = 0.01, alpha = 0.7, size = 2) + 
    geom_text_repel(data = temp_data,
                    aes(label = Variant_name),
                    size = 3,
                    box.padding = 0.5,
                    point.padding = 0.3,
                    max.overlaps = 30) +
    scale_y_log10() + 
    scale_color_manual(values = color_values) +
    labs(x = "Run length (µm)", 
         y = "Dwell time (s)", 
         color = "EEG seizure lifetime") + 
    theme(legend.position = c(0.85, 0.9));

p2 <- ggplot(KIF1A_merge, aes(x = Run_length, color = EEG_seizure_lifetime)) +
    geom_density(linewidth = 1) +
    scale_color_manual(values = color_values) +
    theme_classic() +
    labs(y = "Density") +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_blank());

p3 <- ggplot(KIF1A_merge, aes(x = Dwell_time, color = EEG_seizure_lifetime)) +
    geom_density(linewidth = 1) +
    scale_x_log10() +
    scale_color_manual(values = color_values) +
    coord_flip() +
    theme_classic() +
    labs(y = "Density") +
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          axis.text.y = element_blank());

p4 <- ggplot() + 
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA));

p_merge <- plot_grid(p2, p4, p1, p3, ncol = 2, nrow = 2,
                     rel_widths = c(4, 1), rel_heights = c(1, 4), greedy = TRUE);

filename <- paste(scatter_plot_path, "/Run_length_vs_Dwell_time.png", sep="");
ggsave(filename, p_merge, width=7, height=7, units="in");

# ========= #
# Histogram #
# ========= #

histogram_path <- paste(project_path, "/3_Histogram", sep="");
if (!dir.exists(histogram_path)) {
    dir.create(histogram_path);
}

KIF1A_temp <- KIF1A_merge[KIF1A_merge$Instance == 0 & !is.na(KIF1A_merge$VABS_ABC), ];

p <- wrap_plots(lapply(
    list(
        list(label = "Age < 4",  data = KIF1A_temp[KIF1A_temp$Age_VABS < 4,  ]),
        list(label = "Age < 7",  data = KIF1A_temp[KIF1A_temp$Age_VABS < 7,  ]),
        list(label = "Age > 7",  data = KIF1A_temp[KIF1A_temp$Age_VABS > 7,  ]),
        list(label = "Age > 10", data = KIF1A_temp[KIF1A_temp$Age_VABS > 10, ])
    ),
    function(g) {
        ggplot(g$data, aes(x = VABS_ABC)) +
            geom_histogram(binwidth = 5, fill = "#4393C3", color = "white", alpha = 0.8) +
            xlim(min(KIF1A_temp$VABS_ABC, na.rm = TRUE),
                 max(KIF1A_temp$VABS_ABC, na.rm = TRUE)) +
            labs(title = g$label, x = "VABS ABC score", y = "Count") +
            theme_classic(base_size = 11) +
            theme(plot.title = element_text(size = 11, face = "bold"));
    }
), ncol = 1);

filename <- paste(histogram_path, "/Histogram_VABS_ABC_Final.png", sep="");
ggsave(filename, p, width=8, height=8, units="in");

# === Histogram of different pathogenicity score ===
histogram_patho_score <- function(data, var) {
    ggplot(data, aes(x = get(var))) + theme_classic() + 
        geom_histogram(aes(y = ..density..), color = "lightblue", fill = "lightblue") + 
        geom_density(color = "red", linetype=2) + 
        xlab(get_var_description(KIF1A_merge, var)) + ylab("Density") + 
        theme(plot.title = element_text(hjust = 0.5), 
              axis.text = element_text(color = "black"), 
              axis.ticks = element_line(color = "black"));
}

var_list <- c("CADD_phred", "REVEL", "gMVP", "ESM", "AlphaMissense", "MisFit_D", "MisFit_S");

filename <- paste(histogram_path, "/Histogram_of_pathogenicity_score.png", sep="");
p_list <- lapply(var_list, histogram_patho_score, data = KIF1A_merge);
p_code <- paste("ggarrange(", paste(paste("p_list[[", 1:length(p_list), "]]", sep=""), collapse=", "), ", ncol=2, nrow=4)", sep="");
p_combine <- eval(parse(text=p_code));
ggsave(filename, p_combine, width=6, height=6, units="in");

# === Histogram of MisFit-S in different PCA clusters ===
score_name <- "MisFit_S";
filename <- paste(histogram_path, "/Histogram_of_", score_name, ".PCA_cluster.png", sep="");
p_list <- list();
p_list[[1]] <- histogram_patho_score(KIF1A_merge[!is.na(KIF1A_merge$PCA_cluster) & KIF1A_merge$PCA_cluster == 1, ], score_name) + ggtitle("PCA Cluster 1");
p_list[[2]] <- histogram_patho_score(KIF1A_merge[!is.na(KIF1A_merge$PCA_cluster) & KIF1A_merge$PCA_cluster == 2, ], score_name) + ggtitle("PCA Cluster 2");
p_list[[3]] <- histogram_patho_score(KIF1A_merge[!is.na(KIF1A_merge$PCA_cluster) & KIF1A_merge$PCA_cluster == 3, ], score_name) + ggtitle("PCA Cluster 3");
p_list[[4]] <- histogram_patho_score(KIF1A_merge[!is.na(KIF1A_merge$PCA_cluster) & KIF1A_merge$PCA_cluster == 4, ], score_name) + ggtitle("PCA Cluster 4");
p_list[[5]] <- histogram_patho_score(KIF1A_merge[is.na(KIF1A_merge$PCA_cluster), ], score_name) + ggtitle("Other Variants");
p_code <- paste("ggarrange(", paste(paste("p_list[[", 1:length(p_list), "]]", sep=""), collapse=", "), ", ncol=3, nrow=2)", sep="");
p_combine <- eval(parse(text=p_code));
ggsave(filename, p_combine, width=9, height=5, units="in");

# =========== #
# Violin plot #
# =========== #

vioplot_path <- paste(project_path, "/5_Violin_plot", sep="");
if (!dir.exists(vioplot_path)) {
    dir.create(vioplot_path);
}

# ----------------------------------------------------- #
# Correlation of other variables and EEG/seizure status #
# ----------------------------------------------------- #

vioplot_temp <- function(data, var) {
    min_val <- min(data[, var], na.rm=TRUE);
    max_val <- max(data[, var], na.rm=TRUE);
    lab_pos <- min_val + (max_val - min_val) * 1.1;
    ggplot(data[!is.na(data[, var]) & !is.na(data[, group_name]), ], 
           aes(x = get(group_name), y = get(var), fill = get(group_name))) + 
           geom_violin(trim = FALSE) + geom_boxplot(width = 0.2, fill = "white") + 
           ylim(c(min_val, min_val + (max_val - min_val) * 1.2)) + theme_classic() + 
           xlab(get_var_description(KIF1A_anno, group_name)) + ylab(var) + 
           scale_fill_manual(values = c("#A5A5A5", "#5B9BD5")) + 
           stat_compare_means(method = "t.test", color = "#FF0000", 
                              label.x = 1.3, label.y = lab_pos) + 
           theme(legend.position = "none");
}

# === Correlation of pathogenicity score and EEG/seizure status ===
var_group <- c("Seizure_any_lifetime", "EEG_result_lifetime", "EEG_seizure_lifetime");
var_list <- c("CADD_phred", "REVEL", "gMVP", "ESM", "AlphaMissense", "MisFit_S", "MisFit_D");

for (group_name in var_group) {
    filename <- paste(vioplot_path, "/", group_name, "_vs_pathogenicity_score.png", sep="");
    p_list <- lapply(var_list, vioplot_temp, data = KIF1A_merge);
    p_code <- paste("ggarrange(", paste(paste("p_list[[", 1:length(p_list), "]]", sep=""), collapse=", "), ", ncol=3, nrow=3)", sep="");
    p_combine <- eval(parse(text=p_code));
    ggsave(filename, p_combine, width=9, height=7.5, units="in");
}

# === Correlation of molecular data and EEG/seizure status ===
var_group <- c("Seizure_any_lifetime", "EEG_result_lifetime", "EEG_seizure_lifetime");
var_list <- c("Velocity", "Run_length", "Dwell_time", "Conservation", "Hydro_diff", "pI_diff", "MW_diff");

for (group_name in var_group) {
    filename <- paste(vioplot_path, "/", group_name, "_vs_molecular_data.png", sep="");
    p_list <- lapply(var_list, vioplot_temp, data = KIF1A_merge);
    p_code <- paste("ggarrange(", paste(paste("p_list[[", 1:length(p_list), "]]", sep=""), collapse=", "), ", ncol=3, nrow=3)", sep="");
    p_combine <- eval(parse(text=p_code));
    ggsave(filename, p_combine, width=9, height=7.5, units="in");
}

# === Correlation of VABS score and EEG/seizure status ===
vioplot_temp <- function(data, var, group) {
    min_val <- min(data[, var], na.rm=TRUE);
    max_val <- max(data[, var], na.rm=TRUE);
    lab_pos <- min_val + (max_val - min_val) * 1.1;
    ggplot(data[!is.na(data[, var]) & !is.na(data[, group]), ], 
           aes(x = get(group), y = get(var), fill = get(group))) + 
           geom_violin(trim = FALSE) + geom_boxplot(width = 0.2, fill = "white") + 
           ylim(c(min_val, min_val + (max_val - min_val) * 1.2)) + theme_classic() + 
           xlab(get_var_description(KIF1A_anno, group)) + ylab(var) + 
           scale_fill_manual(values = c("#A5A5A5", "#5B9BD5")) + 
           stat_compare_means(method = "t.test", color = "#FF0000", 
                              label.x = 1.3, label.y = lab_pos) + 
           theme(legend.position = "none");
}

var_group <- c("Seizure_any_lifetime", "EEG_result_lifetime", "EEG_seizure_lifetime");
var_list <- c("VABS_ABC_Final", "VABS_ABC_Diff");

p_list <- list();
for (i in c(1:length(var_group))) {
    for (j in c(1:length(var_list))) {
        index <- i + (j - 1) * 3;
        p_list[[index]] <- vioplot_temp(KIF1A_merge, var_list[j], var_group[i]);
    }
}
filename <- paste(vioplot_path, "/EEG_seizure_vs_VABS.png", sep="");
p_code <- paste("ggarrange(", paste(paste("p_list[[", 1:length(p_list), "]]", sep=""), collapse=", "), ", ncol=3, nrow=2)", sep="");
p_combine <- eval(parse(text=p_code));
ggsave(filename, p_combine, width=9, height=5, units="in");

# Violin plot of VABS scores by PCA_cluster
group_name <- "PCA_cluster";
var_list <- c(
    "VABS_ABC_1", "VABS_Communication_1", "VABS_DLS_1", "VABS_Motor_1", "VABS_Socialization_1", 
    "VABS_ABC_2", "VABS_Communication_2", "VABS_DLS_2", "VABS_Motor_2", "VABS_Socialization_2", 
    "VABS_ABC_Diff", "VABS_Communication_Diff", "VABS_DLS_Diff", "VABS_Motor_Diff", "VABS_Socialization_Diff"
);

vioplot_temp <- function(data, var) {
    temp_data <- data[!is.na(data[, group_name]) & !is.na(data[, var]), ];
    temp_data[, group_name] <- as.factor(temp_data[, group_name]);
    ggplot(temp_data, aes(x=get(group_name), y=get(var), fill=get(group_name))) + 
           theme_classic() + theme(legend.position="none") + 
           xlab("PCA Cluster") + ylab(var) + 
           geom_violin(trim=FALSE) + geom_boxplot(width=0.2, fill="white") + 
           theme(axis.text.x = element_text(angle = 45, hjust = 1));
}

filename <- paste(vioplot_path, "/VABS_vs_PCA_cluster.png", sep="");
p_list <- lapply(var_list, vioplot_temp, data = KIF1A_merge);
p_code <- paste("ggarrange(", paste(paste("p_list[[", 1:length(p_list), "]]", sep=""), collapse=", "), ", ncol=5, nrow=3)", sep="");
p_combine <- eval(parse(text=p_code));
ggsave(filename, p_combine, width=13.3, height=7.5, units="in");

# =================== #
# Random Forest model #
# =================== #

prediction_path <- paste(project_path, "/6_Prediction", sep="");
if (!dir.exists(prediction_path)) {
    dir.create(prediction_path);
}

KIF1A_temp <- KIF1A_merge.map[KIF1A_merge.map$Instance == 0, ];
KIF1A_temp[is.na(KIF1A_temp$Seizure_any), ]$Seizure_any <- 0;

temp_data <- KIF1A_merge[KIF1A_merge$Instance == 1, c("Record_ID", "Age_VABS", "VABS_ABC")];
colnames(temp_data) <- c("Record_ID", "Age_VABS_1st");

KIF1A_temp <- merge(KIF1A_temp, temp_data, by = "Record_ID", all.x = TRUE, all.y = FALSE);
KIF1A_temp <- KIF1A_temp[!(KIF1A_temp$Record_ID %in% removed_IDs), ];
KIF1A_temp$GBTM_cluster <- factor(KIF1A_temp$GBTM_cluster); # Use factor for classification

var_list <- c("Optic_nerve_change", "Hypertonia", "Hypotonia", "Microcephaly", "Ataxia", 
              "Dystonia", "Neuropathy", "Developmental_concern", "gMVP", "AlphaMissense", "MisFit_S", 
              "Conservation", "Hydro_diff", "pI_diff", "MW_diff");

# ----------------------------------------------------------- #
# Model 1: VABS ABC ~ Age + clin_vars + patho_vars + mol_vars #
# ----------------------------------------------------------- #

predictor_list <- list(
    "Full"    = c("Age_VABS", clin_vars, patho_vars, mol_vars), # Full model
    "R1_age"  = c(clin_vars, patho_vars, mol_vars),             # Ablation (Remove age)
    "R2_mol"  = c("Age_VABS", clin_vars, patho_vars),           # Ablation (Remove molecular)
    "R3_both" = c(clin_vars, patho_vars)                        # Ablation (Remove both)
);
response <- "VABS_ABC";

# === Random Forest model ===
VABS_ABC_Final.RF <- list();
filename <- file.path(prediction_path, "VABS_ABC_Final.RF.rds");
if (!file.exists(filename)) {
    for (name in names(predictor_list)) {
        message("Training ", response, " model: ", name);
        predictor <- predictor_list[[name]];
        filepath <- file.path(prediction_path, paste("VABS_ABC_Final.RF.", name, sep=""));
        model <- regression_model(
            data      = KIF1A_temp[!is.na(KIF1A_temp$VABS_ABC) & 
                                   !is.na(KIF1A_temp$Molecular_cluster), ], 
            predictor = predictor, 
            response  = response, 
            model     = "rf", 
            id_col    = "Record_ID", 
            filepath  = filepath
        );
        VABS_ABC_Final.RF[[name]] <- model;
    }
    saveRDS(VABS_ABC_Final.RF, filename);
} else {
    VABS_ABC_Final.RF <- readRDS(filename);
}

filepath <- file.path(prediction_path, "VABS_ABC_Final.RF.Evaluation");
model_evaluation(VABS_ABC_Final.RF, KIF1A_temp, "Age_VABS", 10, 20, "Final VABS ABC", filepath);

# -------------------------------------------------------------------------- #
# Model 2: GBTM_cluster ~ Age_1st + ΔAge + clin_vars + patho_vars + mol_vars #
# -------------------------------------------------------------------------- #

predictor_list <- list(
    "Full"    = c("Age_VABS", clin_vars, patho_vars, mol_vars), # Full model
    "R1_age"  = c(clin_vars, patho_vars, mol_vars),             # Ablation (Remove age)
    "R2_mol"  = c("Age_VABS", clin_vars, patho_vars),           # Ablation (Remove molecular)
    "R3_both" = c(clin_vars, patho_vars)                        # Ablation (Remove both)
);
response <- "GBTM_cluster";

# === Random Forest model (classification) ===
GBTM_cluster.RF <- list();
filename <- file.path(prediction_path, "GBTM_cluster.RF.rds");
if (!file.exists(filename)) {
    for (name in names(predictor_list)) {
        message("Training ", response, " model: ", name);
        predictor <- predictor_list[[name]];
        filepath  <- file.path(prediction_path, paste("GBTM_cluster.RF.", name, sep = ""));
        model <- classification_model(
            data        = KIF1A_temp[!is.na(KIF1A_temp$GBTM_cluster) &
                                     !is.na(KIF1A_temp$Molecular_cluster), ],
            predictor   = predictor,
            response    = response,
            model       = "rf",
            id_col      = "Record_ID",
            filepath    = filepath
        );
        GBTM_cluster.RF[[name]] <- model;
    }
    saveRDS(GBTM_cluster.RF, filename);
} else {
    GBTM_cluster.RF <- readRDS(filename);
}

filepath <- file.path(prediction_path, "GBTM_cluster.RF.Evaluation");
model_evaluation(GBTM_cluster.RF, KIF1A_temp, type = "classification", filepath = filepath);

# ------------------------------------------------------------------- #
# Model 3: ΔVABS ~ Age_1st + ΔAge + clin_vars + patho_vars + mol_vars #
# ------------------------------------------------------------------- #

predictor_list <- list(
    "Full"    = c("Age_VABS_1st", "Age_VABS_Diff", clin_vars, patho_vars, mol_vars), # Full model
    "R1_age"  = c("Age_VABS_Diff", clin_vars, patho_vars, mol_vars),                 # Ablation (Remove age_1st)
    "R2_mol"  = c("Age_VABS_1st", "Age_VABS_Diff", clin_vars, patho_vars),           # Ablation (Remove molecular)
    "R3_both" = c("Age_VABS_Diff", clin_vars, patho_vars)                            # Ablation (Remove both)
);
response <- "VABS_ABC_Diff";

# === Random Forest model ===
VABS_ABC_Diff.RF <- list();
filename <- file.path(prediction_path, "VABS_ABC_Diff.RF.rds");
if (!file.exists(filename)) {
    for (name in names(predictor_list)) {
        message("Training ", response, " model: ", name);
        predictor <- predictor_list[[name]];
        filepath <- file.path(prediction_path, paste("VABS_ABC_Diff.RF.", name, sep=""));
        model <- regression_model(
            data      = KIF1A_temp[!is.na(KIF1A_temp$Age_VABS_Diff) & 
                                   !is.na(KIF1A_temp$Molecular_cluster), ], 
            predictor = predictor, 
            response  = response, 
            model     = "rf", 
            id_col    = "Record_ID", 
            filepath  = filepath
        );
        VABS_ABC_Diff.RF[[name]] <- model;
    }
    saveRDS(VABS_ABC_Diff.RF, filename);
} else {
    VABS_ABC_Diff.RF <- readRDS(filename);
}

filepath <- file.path(prediction_path, "VABS_ABC_Diff.RF.Evaluation");
model_evaluation(VABS_ABC_Diff.RF, KIF1A_temp, NULL, NULL, 5, "Delta VABS ABC", filepath);

# ------------------------------------------------------------------ #
# Model 4: Rate ~ Age_1st + ΔAge + clin_vars + patho_vars + mol_vars #
# ------------------------------------------------------------------ #

predictor_list <- list(
    "Full"    = c("Age_VABS_1st", "Age_VABS_Diff", clin_vars, patho_vars, mol_vars), # Full model
    "R1_age"  = c("Age_VABS_Diff", clin_vars, patho_vars, mol_vars),                 # Ablation (Remove age_1st)
    "R2_mol"  = c("Age_VABS_1st", "Age_VABS_Diff", clin_vars, patho_vars),           # Ablation (Remove molecular)
    "R3_both" = c("Age_VABS_Diff", clin_vars, patho_vars)                            # Ablation (Remove both)
);
response <- "VABS_ABC_Rate";

# === Random Forest model ===
VABS_ABC_Rate.RF <- list();
filename <- file.path(prediction_path, "VABS_ABC_Rate.RF.rds");
if (!file.exists(filename)) {
    for (name in names(predictor_list)) {
        message("Training ", response, " model: ", name);
        predictor <- predictor_list[[name]];
        filepath <- file.path(prediction_path, paste("VABS_ABC_Rate.RF.", name, sep=""));
        model <- regression_model(
            data      = KIF1A_temp[!is.na(KIF1A_temp$Age_VABS_Diff) & 
                                   !is.na(KIF1A_temp$Molecular_cluster), ], 
            predictor = predictor, 
            response  = response, 
            model     = "rf", 
            id_col    = "Record_ID", 
            filepath  = filepath
        );
        VABS_ABC_Rate.RF[[name]] <- model;
    }
    saveRDS(VABS_ABC_Rate.RF, filename);
} else {
    VABS_ABC_Rate.RF <- readRDS(filename);
}

filepath <- file.path(prediction_path, "VABS_ABC_Rate.RF.Evaluation");
model_evaluation(VABS_ABC_Rate.RF, KIF1A_temp, NULL, NULL, 5, "Rate of Change", filepath);









# Elastic Net model
filename <- paste(prediction_path, "/VABS_ABC_Final.EN.rds", sep="");
if (!file.exists(filename)) {
    VABS_ABC_Final.EN <- regression_model(
        data         = KIF1A_temp, 
        predictor    = predictor,
        response     = response,
        model        = "en",
        id_col       = "Record_ID",
        filepath     = paste(prediction_path, "/VABS_ABC_Final.EN", sep=""),
        tune_folds   = 5,
        tune_repeats = 3,
        n_cores      = 20
    );
    saveRDS(VABS_ABC_Final.EN, filename);
} else {
    VABS_ABC_Final.EN <- readRDS(filename);
}

# XGBoost model
filename <- paste(prediction_path, "/VABS_ABC_Final.XGB.rds", sep="");
if (!file.exists(filename)) {
    VABS_ABC_Final.XGB <- regression_model(
        data      = KIF1A_temp,
        predictor = predictor,
        response  = response,
        model     = "xgb",
        id_col    = "Record_ID",
        filepath  = paste(prediction_path, "/VABS_ABC_Final.XGB", sep=""),
        n_cores   = 20
    );
    saveRDS(VABS_ABC_Final.XGB, filename);
} else {
    VABS_ABC_Final.XGB <- readRDS(filename);
}


# ---------------------------------------------------------------------------- #
# Model: VABS ~ Sex + Age + EEG/seizure + Pathogenicity score + Molecular data #
# ---------------------------------------------------------------------------- #
predictor <- c("Sex", "Age_VABS", "Optic_nerve_change", "Seizure_any", "EEG_result", "Hypertonia", "Hypotonia", 
    "Microcephaly", "Ataxia", "Dystonia", "Neuropathy", "REVEL", "ESM", "MisFit_D", "Movement", "Developmental_concern", 
    "Diffusion", "Velocity", "Run_length", "Dwell_time", "Conservation", "Hydro_diff", "pI_diff", "MW_diff");
response <- "VABS_ABC";

filepath <- paste(prediction_path, "/VABS_ABC_Final.new", sep="");
RF_model.VABS_ABC_Final <- RF_model(KIF1A_temp, predictor, response, 
    model_name = "RF_model.top_freq", model_type = "regression", 
    split_method = "top_freq", top_count = 10, filepath = filepath);

# Plot model prediction
temp_data <- RF_model.VABS_ABC_Final$prediction;
temp_data$Residual <- abs(temp_data$Actual - temp_data$Predicted);
temp_data$Variant_name <- gsub("^p.", "", temp_data$P_notation);
temp_data[temp_data$Residual <= 10, ]$Variant_name <- "";
r2_value <- cor(temp_data$Actual, temp_data$Predicted) ^ 2;

p <- ggplot(temp_data, aes(x = Actual, y = Predicted)) + 
    theme_classic() + xlim(20, 140) + ylim(20, 140) + 
    geom_point(aes(color = Residual), size = 4, alpha = 0.8) +
    scale_color_gradient(low = "#1A9641", high = "#D7191C",
                         name = "Absolute Error") + 
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    geom_text_repel(aes(label = Variant_name),
                    size = 3,
                    box.padding = 0.5,
                    max.overlaps = 30) + 
    annotate("text", x = 25, y = 90, 
             label = paste("R^2 == ", sprintf("%0.2f", r2_value), sep = ""), 
             parse = TRUE, 
             size = 5, fontface = "bold") +
    labs(x = "Actual Final VABS ABC", 
         y = "Predicted Final VABS ABC") + 
    theme(legend.position = c(0.90, 0.15));
filename <- paste(filepath, "/RF_model.top_freq.prediction.png", sep="");
ggsave(filename, p, width=6, height=6, units="in");

# --------------------------------------------------------------------- #
# Model: EEG/seizure ~ Sex + Age + Pathogenicity score + Molecular data #
# --------------------------------------------------------------------- #
predictor <- c("Sex", "Age_VABS", "REVEL", "ESM", "MisFit_D", "Movement", "Diffusion", "Velocity", 
    "Run_length", "Dwell_time", "Conservation", "Hydro_diff", "pI_diff", "MW_diff");
response <- "EEG_or_seizure";

KIF1A_temp <- KIF1A_merge.map[KIF1A_merge.map$Instance == 0, ];
KIF1A_temp[KIF1A_temp$EEG_or_seizure == 1, ]$EEG_or_seizure <- 0;
KIF1A_temp[KIF1A_temp$EEG_or_seizure == 2, ]$EEG_or_seizure <- 1;

filepath <- paste(prediction_path, "/EEG_or_seizure", sep="");
RF_model.EEG_seizure_lifetime <- RF_model(KIF1A_temp, predictor, response, 
    model_name = "RF_model.top_freq", model_type = "classification", 
    split_method = "top_freq", top_count = 10, filepath = filepath);


filepath <- paste(prediction_path, "/EEG_or_seizure.age_above_7", sep="");
KIF1A_temp <- KIF1A_merge.map[!is.na(KIF1A_merge.map$Age_VABS_Final) & KIF1A_merge.map$Age_VABS_Final > 7, ];
RF_model.EEG_seizure_lifetime <- RF_model(KIF1A_temp, predictor, response, 
    model_name = "RF_model.top_freq", model_type = "classification", 
    split_method = "top_freq", top_count = 10, filepath = filepath);

temp_data <- RF_model.EEG_seizure_lifetime$prediction;
actual_binary <- ifelse(temp_data$Actual == "Class_1", 1, 0)

roc_obj <- roc(actual_binary, temp_data$Class_1)
auc_value <- auc(roc_obj)

roc_data <- data.frame(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities
)

p <- ggplot(roc_data, aes(x = FPR, y = TPR)) +
    geom_line(color = "#2E86AB", size = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    annotate("text", x = 0.6, y = 0.2, 
             label = paste0("AUC = ", round(auc_value, 3)), 
             size = 5, fontface = "bold") +
    labs(x = "False Positive Rate", 
         y = "True Positive Rate",
         title = "ROC Curve for EEG Seizure Prediction") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(paste0(filepath, "/", "RF_model.top_freq.ROC.png"), 
       plot = p, width = 6, height = 6, dpi = 300)

# ====================== #
# Clinical data analysis #
# ====================== #

clinical_results_path <- paste(project_path, "/8_Clinical_results", sep="");
if (!dir.exists(clinical_results_path)) {
    dir.create(clinical_results_path);
}

KIF1A_temp <- KIF1A_merge.map[KIF1A_merge.map$Instance == 0, ];

temp_data <- KIF1A_merge[KIF1A_merge$Instance == 1, c("Record_ID", "Age_VABS")];
colnames(temp_data) <- c("Record_ID", "Age_VABS_1st");

KIF1A_temp <- merge(KIF1A_temp, temp_data, by = "Record_ID", all.x = TRUE, all.y = FALSE);
KIF1A_temp <- KIF1A_temp[!(KIF1A_temp$Record_ID %in% removed_IDs), ];

# ---------------------------------- #
# Molecular phenotype vs EEG/seizure #
# ---------------------------------- #

out_1 <- logistic_regression(KIF1A_merge.map, KIF1A_anno, D_var = "EEG_seizure_lifetime", I_var = mol_vars, 
                             C_var = NA, group = NA, accuracy = 2, adj_only = FALSE, min_allow = 3);

out_list <- list(out_1);
clinical_analysis_report(project_name="Molecular_phenotype_vs_EEG_seizure_lifetime", out_list=out_list, font_size=15);

# ----------------------------------------------- #
# Baseline clinical features of each GBTM_cluster #
# ----------------------------------------------- #

temp_data <- KIF1A_temp;
var_list <- c("Sex", "Age_VABS", "VABS_ABC", clin_vars);
for (var in setdiff(var_list, c("Sex", "Age_VABS", "VABS_ABC", "Seizure_any"))) {
    temp_data[temp_data[, var] == 1, var] <- 0;
    temp_data[temp_data[, var] == 2, var] <- 1;
}
out_1 <- descriptive_statistics(temp_data, KIF1A_anno, var_list, group = "GBTM_cluster");

# ------------------------------------------------ #
# Molecular clusters and GBTM trajectory subgroups #
# ------------------------------------------------ #
out_2 <- descriptive_statistics(KIF1A_temp, KIF1A_anno, var_list = "Molecular_cluster", group = "GBTM_cluster");

# ------------------------------------------------ #
# Linear regression of multiple variables and VABS #
# ------------------------------------------------ #
I_var <- c(clin_vars, patho_vars, mol_vars);

temp_data <- KIF1A_temp;
for (var in I_var[c(1, 3:10)]) {
    temp_data[temp_data[, var] == 1, var] <- 0;
    temp_data[temp_data[, var] == 2, var] <- 1;
}

# === multiple variables and final VABS ABC (age >= 10) ===
D_var <- "VABS_ABC";
C_var <- "Age_VABS";
out_3 <- linear_regression(temp_data[temp_data$Age_VABS >= 10, ], KIF1A_anno, D_var, I_var, C_var, 
                           group=NA, accuracy=2, adj_only=FALSE, normal_test=FALSE);

# === multiple variables and change of VABS ABC ===
D_var <- "VABS_ABC_Diff";
C_var <- c("Age_VABS_1st", "Age_VABS_Diff");
out_4 <- linear_regression(temp_data, KIF1A_anno, D_var, I_var, C_var, group=NA, accuracy=2, adj_only=TRUE, normal_test=FALSE);

# === multiple variables and rate of change of VABS ABC ===
D_var <- "VABS_ABC_Rate";
C_var <- c("Age_VABS_1st", "Age_VABS_Diff");
out_5 <- linear_regression(temp_data, KIF1A_anno, D_var, I_var, C_var, group=NA, accuracy=2, adj_only=TRUE, normal_test=FALSE);

# ------------------------------------------------------------ #
# Logistic regression of molecular data and seizure resistance #
# ------------------------------------------------------------ #
temp_data <- KIF1A_temp[!is.na(KIF1A_temp$Molecular_cluster) & KIF1A_temp$Seizure_resistant != 1, ];
temp_data[temp_data$Seizure_resistant == 2, "Seizure_resistant"] <- 1;

I_var <- mol_vars;
C_var <- c("Sex", "Age_VABS");
D_var <- "Seizure_resistant";
out_6 <- logistic_regression(temp_data, KIF1A_anno, D_var, I_var, C_var, group=NA, accuracy=2, adj_only=FALSE, min_allow=3);

out_list <- list(out_1, out_2, out_3, out_4, out_5, out_6);
clinical_analysis_report(project_name="Clinical_analysis_results", out_list=out_list, font_size=15);

# ================= #
# Save the Rsession #
# ================= #

filename <- paste(project_path, "/KIF1A_", gsub("-", "", Sys.Date(), fixed=TRUE), ".RData", sep="");
save.image(filename);