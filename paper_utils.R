suppressMessages({
#     source("/data/srlab/ik936/Foxxy/utils/utils.R")
    library(ggthemes)
    library(magrittr)
    library(Matrix)
    library(dplyr)
    library(ggrepel)
    library(tidyr)
    library(data.table)
    
    library(RColorBrewer)
    library(scales)
    library(patchwork)  
})

library(Rcpp)
library(RcppArmadillo)


fig.size <- function (height, width) {
    options(repr.plot.height = height, repr.plot.width = width)
}


## our libraries
library(reticulate)
library(ggrastr)
do_scatter <- function(umap_use, meta_data, label_name, no_guides = TRUE, do_labels = TRUE, nice_names, palette_use,
                       pt_size = 4, point_size = .5, base_size = 12, do_points = TRUE, do_density = FALSE, h = 6, w = 8) {
    plt_df <- umap_use %>% data.frame() %>% 
        cbind(meta_data) %>% 
        dplyr::sample_frac(1L) 
    plt_df$given_name <- plt_df[[label_name]]
    
    if (!missing(nice_names)) {
        plt_df %<>%
            dplyr::inner_join(nice_names, by = "given_name") %>% 
            subset(nice_name != "" & !is.na(nice_name))

        plt_df[[label_name]] <- plt_df$nice_name        
    }
        
    plt <- plt_df %>% 
        ggplot(aes_string("X1", "X2", col = label_name, fill = label_name)) + 
        theme_tufte(base_size = base_size) + 
        theme(panel.background = element_rect(fill = NA, color = "black")) + 
        guides(color = guide_legend(override.aes = list(stroke = 1, alpha = 1, shape = 16, size = 4)), alpha = FALSE) +
        scale_color_manual(values = palette_use) + 
        scale_fill_manual(values = palette_use) +    
        theme(plot.title = element_text(hjust = .5)) + 
        labs(x = "UMAP 1", y = "UMAP 2") 
    
    if (do_points) 
        plt <- plt + geom_point_rast(dpi = 300, width = w, height = h, size = point_size) 
#         plt <- plt + geom_point()
    if (do_density) 
        plt <- plt + geom_density_2d()    
        

    if (no_guides)
        plt <- plt + guides(col = FALSE, fill = FALSE, alpha = FALSE)
    
    if (do_labels) {
        plt <- plt + geom_text_repel(
            data = data.table(plt_df)[, .(X1 = mean(X1), X2 = mean(X2)), by = label_name], 
            label.size = NA, aes_string(label = label_name), 
            color = "black", size = pt_size, alpha = 1, segment.size = 0) + 
        guides(col = FALSE, fill = FALSE)
        
    }
    return(plt)
}

plot_clusters3 <- function(cluster_ids, labels, pt_size = 14, umap_use = umap_post, 
                           do_labels = FALSE, palette_use = tableau_color_pal("Tableau 20")(20), 
                           do_points = TRUE, do_density = FALSE, alpha_pt = 1, min_cluster_size = 20) {
    cluster_table <- table(cluster_ids)
    clusters_keep <- names(which(cluster_table > min_cluster_size))
    plt_df <- umap_use %>% data.frame() %>% cbind(cluster = cluster_ids) %>%
        subset(cluster %in% clusters_keep) 
    plt <- plt_df %>% 
        dplyr::sample_frac(1L) %>%
        ggplot(aes(X1, X2, col = factor(cluster)))
    if (do_points) 
        plt <- plt + geom_point(shape = '.', alpha = alpha_pt)
    if (do_density) 
        plt <- plt + geom_density_2d()    
    
    
    plt <- plt + 
        theme_tufte() + geom_rangeframe(col = "black") + 
#         theme(axis.line = element_line()) +
        guides(color = guide_legend(override.aes = list(stroke = 1, alpha = 1, shape = 16, size = 4))) + 

        scale_color_manual(values = palette_use) + 
        labs(x = "UMAP 1", y = "UMAP 2") +
        theme(plot.title = element_text(hjust = .5))
    
    if (do_labels) 
        plt <- plt + geom_label_repel(data = data.table(plt_df)[, .(X1 = mean(X1), X2 = mean(X2)), by = cluster], 
                                aes(label = cluster), size = pt_size, alpha = .8) + 
        guides(col = FALSE)
    return(plt)
}



setupVals <- function(data_mat, feature, qlo, qhi) {
    .x <- data_mat[feature, , drop = FALSE] %>% as("dgTMatrix")
    cutoffs <- quantileSparse(.x, c(qlo, qhi))
    cutoffs[2] <- max(cutoffs[2], min(.x@x))
    if (qlo == 0 & qhi == 1) {
        return(.x)
    } 
    
    if (qlo > 0) {
        .x@x[.x@x < cutoffs[1]] <- cutoffs[1]
#         message(sprintf("For %s, lo = %.3f", feature, ifelse(length(.x@x) == ncol(.x), cutoffs[1], NA)))
    }
    if (qhi < 1) {
        .x@x[.x@x > cutoffs[2]] <- cutoffs[2]
#         message(sprintf("For %s, hi = %.3f", feature, cutoffs[2]))
        
    }
    return(.x)
}


quantileSparse <- function(.x, qlist) {
    ratio_zero <- 1 - (length(.x@x) / ncol(.x))
    q_nz <- which(qlist > ratio_zero)
    q_adj <- (qlist[q_nz] - ratio_zero) / (1 - ratio_zero)
    res <- rep(0, length(qlist))
    res[q_nz] <- quantile(.x@x, q_adj)
    res
}

## TODO: test is feature is present
## TODO: allow for different cutoffs, for each marker
## TODO: somehow draw canvas first, then do plotting? 
library(patchwork)
library(ggthemes)

plotFeatures <- function(data_mat, dim_df, features, nrow = 1, 
                         qlo = 0.05, qhi = 1, order_by_expression = FALSE, 
                         pt_shape = 16, pt_size = .5, no_guide = FALSE,
                         .xlim = c(NA, NA), .ylim = c(NA, NA), color_high = muted("blue")) {
    plt_df <- data.frame(dim_df[, 1:2])
    colnames(plt_df) <- c("X1", "X2")


    plt_list <- lapply(features, function(feature) {
        .x <- setupVals(data_mat, feature, qlo, qhi)
        plt_df$value <- 0
        plt_df[.x@j + 1, "value"] <- .x@x
        if (order_by_expression) {
            plt_df %<>% dplyr::arrange(value)             
        } else {
            plt_df %<>% dplyr::sample_frac(1L)
        }

        plt <- plt_df %>% 
            ggplot(aes(X1, X2, color = value)) + 
            geom_point_rast(dpi = 300, width = 6, height = 4, size = .5, shape = pt_shape) + 
#             geom_point(shape = ".") + 
            scale_color_gradient2(na.value = "lightgrey", mid = "lightgrey", midpoint = 0, high = color_high) + 
            theme_tufte(base_size = 14, base_family = "Helvetica") + 
            theme(panel.background = element_rect(), plot.title = element_text(hjust = .5)) +
            labs(x = "UMAP 1", y = "UMAP 2", title = feature) + 
            NULL
        if (no_guide) {
            plt <- plt + 
            guides(color = FALSE) 
        }
        
        if (sum(is.na(.xlim)) < 2) 
            plt <- plt + xlim(.xlim)
        if (sum(is.na(.ylim)) < 2) 
            plt <- plt + ylim(.ylim)
        plt

    })
    if (length(plt_list) > 1) {
        Reduce(`+`, plt_list) + patchwork::plot_layout(nrow = nrow)
    } else {
        plt_list[[1]]
    }
}



# fig.size(5, 10)
do_dot_plot <- function(gene_levels) {
    data_plot <- markers_tidy %>%
        subset(symbol %in% gene_levels) %>%
        dplyr::mutate(symbol = factor(symbol, gene_levels)) %>%
        dplyr::mutate(pct = ifelse(is.na(pct), 0, pct)) %>% 
        dplyr::mutate(avg_exp = ifelse(is.na(avg_exp), 0, avg_exp)) %>% 
        dplyr::group_by(symbol) %>% 
        mutate(avg_exp = scale(avg_exp)) %>% 
        mutate(avg_exp_scale = as.numeric(cut(MinMax(avg_exp, max = 2.5, min = -2.5), breaks = 10))) %>%
        ungroup() %>%
        identity()
    data_plot <- data.table(data_plot)[, ptcolor := colorRampPalette(c("grey", "blue"))(10)[avg_exp_scale], by = .(cluster, symbol)]    
    
    data_plot %>% 
    #     dplyr::mutate(cell_type = factor(cell_type, rev(group_levels))) %>% 
        ggplot(aes(reorder(cluster, as.integer(gsub("X", "", cluster))), symbol)) + 
#         ggplot(aes(cluster, symbol)) + 
            geom_point(aes(size = pct, color = ptcolor)) + 
            scale_size_area(max_size = 4) + 
    #         scale_radius(range = c(0, 4)) + 
            scale_color_identity() + 
            theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
            theme_tufte() + theme(panel.background = element_rect()) + 
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
            labs(x = "", y = "") + 
            NULL
    
}

MinMax <- function (data, min, max) {
    data2 <- data
    data2[data2 > max] <- max
    data2[data2 < min] <- min
    return(data2)
}

sourceCpp("/data/srlab/ik936/fast_wilcox.cpp")
source("/data/srlab/ik936/Foxxy/utils/fast_wilcox.R")

## for sumGroups
sumOverRowNames <- function(X) {
    name_factors <- factor(row.names(X))
    res <- sumGroups(X@x, X@p, X@i, ncol(X), as.integer(name_factors) - 1, length(levels(name_factors)))
    row.names(res) <- levels(name_factors)[1:nrow(res)]
    colnames(res) <- colnames(X)
    return(res)
}

read10x <- function(run, suffix) {
    barcode.loc <- file.path(run, "barcodes.tsv")
    gene.loc <- file.path(run, "genes.tsv")
    matrix.loc <- file.path(run, "matrix.mtx")

    data <- readMM(file = matrix.loc) %>% as("dgCMatrix")
    cell.names <- readLines(barcode.loc)
    cell.names <- gsub("-1$", "", cell.names)
    if (!missing(suffix)) {
        cell.names %<>% paste(suffix, sep = "_")
    }
    
    gene.names <- fread(gene.loc, header = FALSE)$V2
    row.names(data) <- gene.names
    colnames(data) <- cell.names

    
    return(as(sumOverRowNames(data), "dgCMatrix"))
}
    
cosine_normalize <- function(X, MARGIN = 1, do_safe = TRUE) {
    if (do_safe) {
        X <- sweep(X, MARGIN, apply(X, MARGIN, max), "/")
    }
    sweep(X, MARGIN, apply(X, MARGIN, function(x) sqrt(sum(x^2))), "/")
}

column_to_rownames <- function(df, colname) {
    df %<>% data.frame()
    row.names(df) <- df[[colname]]
    df[[colname]] <- NULL
    return(df)
}                           
