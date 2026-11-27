#!/usr/bin/env Rscript
# Publication figures using base R and Cairo. Run Rscript code/output/plot_analysis.R.

script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
root <- dirname(normalizePath(script, winslash="/"))
requested_figures <- commandArgs(trailingOnly=TRUE)
figures <- file.path(root, "figures")
dir.create(figures, showWarnings=FALSE)
table <- function(name) read.delim(file.path(root, "tables", paste0(name, ".tsv")),
                                 na.strings=c("", "nan"), check.names=FALSE)
families <- read.delim(file.path(root, "..", "rvne", "results", "supplemental_run", "families.tsv"),
                       na.strings=c("", "nan"), check.names=FALSE)
levels <- table("table_by_setting")
profiles <- table("figure_position_profiles")
curves <- table("figure_d_curves")
gap <- table("table_reference_gap")

ink <- "#152333"
muted <- "#516173"
gridline <- "#E2E8EE"
colors <- c(E0="#586677", E2="#CE701B", E3="#9E5615", E4="#8793A2",
            E5="#744B98", E6="#926D53", E7="#007F6D", E8="#1D57A0",
            Shannon="#7D8997", "AL2CO variance"="#51687B")
labels <- c(E0="E0 Original", E4="E4 Adaptive", E5="E5 Segmentation",
            E6="E6 Default SCAW", E7="E7 Retuned RVNE", E8="E8 Retuned SCAW",
            Shannon="Shannon", "AL2CO variance"="AL2CO variance")

save_figure <- function(name, width, height, draw) {
    if (length(requested_figures) && !name %in% requested_figures) return(invisible(NULL))
    render <- function() {
        on.exit(dev.off())
        par(family="sans", fg=ink, col.axis=ink, col.lab=ink, col.main=ink, las=1,
            tcl=-0.25, mgp=c(2.05, 0.55, 0), lend="round", ljoin="round")
        draw()
    }
    png(file.path(figures, paste0(name, ".png")), width=width, height=height,
        units="in", res=300, type="cairo", bg="white", pointsize=11)
    render()
    cairo_pdf(file.path(figures, paste0(name, ".pdf")), width=width, height=height,
              bg="white", pointsize=11)
    render()
}

panel <- function(letter, title) {
    mtext(letter, side=3, line=1.17, adj=0, font=2, cex=1.14, col=ink)
    mtext(title, side=3, line=1.20, adj=0, at=par("usr")[1] +
          0.085 * diff(par("usr")[1:2]), font=2, cex=0.99, col=ink)
}

axis_frame <- function(x, y, xlabels=x, ylabels=y, grid="y") {
    if (grid == "x") abline(v=x, col=gridline, lwd=0.85)
    if (grid == "y") abline(h=y, col=gridline, lwd=0.85)
    axis(1, at=x, labels=xlabels, col=ink, cex.axis=0.91)
    axis(2, at=y, labels=ylabels, col=ink, cex.axis=0.91)
    box(bty="l", col=ink)
}

values_for <- function(pattern, method, metric) {
    families[families$scenario == pattern & families$method == method, metric]
}

draw_distribution <- function(values, at, color, horizontal=FALSE) {
    boxplot(values, at=at, horizontal=horizontal, add=TRUE, axes=FALSE, outline=FALSE,
            boxwex=0.46, col=adjustcolor(color, alpha.f=0.16), border=color,
            pars=list(medlwd=1.6, staplelty=1, whisklty=1))
    offset <- seq(-0.155, 0.155, length.out=length(values))
    if (horizontal) {
        points(values, at+offset, pch=16, col=adjustcolor(color, alpha.f=0.56), cex=0.69)
    } else {
        points(at+offset, values, pch=16, col=adjustcolor(color, alpha.f=0.56), cex=0.69)
    }
}

save_figure("fig01_patterns_abc", 12.4, 3.95, function() {
    layout(matrix(1:3, nrow=1), widths=c(1, 1.19, 1.27))
    par(mar=c(3.5, 4.2, 3.1, 1.05))
    methods <- c("E5", "E6", "E8")
    plot(NA, xlim=c(0.52, 3.63), ylim=c(-0.19, 4.36), xaxs="i", yaxs="i",
         xlab="", ylab="False boundaries per family", axes=FALSE)
    axis_frame(1:3, 0:4, methods)
    for (i in seq_along(methods)) {
        method <- methods[i]
        color <- colors[method]
        values <- values_for("A", method, "false_boundaries")
        value <- mean(values)
        rect(i-0.23, 0, i+0.23, value, col=adjustcolor(color, alpha.f=0.76), border=NA)
        offset <- seq(-0.155, 0.155, length.out=length(values))
        points(i+offset, values, pch=21, bg=adjustcolor(color, alpha.f=0.64),
               col="white", lwd=0.55, cex=0.72)
        text(i, value+0.07, sprintf("%.2f", value), pos=3, offset=0.12, cex=0.85, col=ink)
    }
    panel("A", "Uniform conservation")

    par(mar=c(3.5, 4.2, 3.1, 1.0))
    methods <- c("E0", "E2", "E3", "E5", "E6", "E7", "E8")
    plot(NA, xlim=c(0.5, 7.5), ylim=c(0.93, 1.008), xaxs="i", yaxs="i",
         xlab="", ylab="Variable-site IoU", axes=FALSE)
    ticks <- seq(0.94, 1.00, 0.02)
    axis_frame(seq_along(methods), ticks, methods, sprintf("%.2f", ticks))
    for (i in seq_along(methods)) {
        method <- methods[i]
        draw_distribution(values_for("B", method, "variable_site_iou"), i, colors[method])
    }
    panel("B", "Central variable block")

    par(mar=c(3.5, 4.2, 3.1, 0.85))
    plot(NA, xlim=c(0.55, 3.45), ylim=c(-0.06, 1.22), xlab="Motif width (residues)",
         ylab="Conserved motif recall", axes=FALSE, xaxs="i", yaxs="i")
    axis_frame(1:3, seq(0,1,0.25), c(3,5,9), paste0(seq(0,100,25), "%"))
    methods <- c("E2", "E3", "E8")
    for (i in seq_along(methods)) {
        method <- methods[i]
        data <- levels[levels$pattern == "C" & levels$method == method & levels$metric == "sensitivity", ]
        x <- match(data$level, c(3,5,9)) + c(-0.24,0,0.24)[i]
        rect(x-0.095, 0, x+0.095, data$mean, col=colors[method], border=NA)
        nonzero <- data$mean > 0
        text(x, data$mean+0.038, sprintf("%g%%", round(100*data$mean, 1)),
             cex=0.70, col=ifelse(nonzero, ink, muted))
    }
    legend("top", c("E2 Uniform", "E3 Gaussian", "E8 Retuned"), bty="n", cex=0.79,
           fill=colors[methods], border=NA, horiz=TRUE,
           x.intersp=0.55, text.col=ink)
    panel("C", "Short conserved motif")
})

d_panel <- function(methods, metric, limits, ticks, tick_labels, xlab, letter, title, stem=FALSE) {
    y <- rev(seq_along(methods))
    plot(NA, xlim=limits, ylim=c(0.45, length(y)+1.05), xlab=xlab, ylab="",
         axes=FALSE, xaxs="i", yaxs="i")
    segments(ticks, 0.45, ticks, length(y)+0.39, col=gridline, lwd=0.85)
    labels_x <- limits[2] - 0.013 * diff(limits)
    text(labels_x, length(y)+0.67, "Mean", adj=1, cex=0.72, font=2, col=muted)
    for (i in seq_along(methods)) {
        method <- methods[i]
        color <- colors[method]
        values <- values_for("D", method, metric)
        average <- mean(values)
        hero <- method == if (stem) "E8" else "E5"
        if (stem) {
            segments(0, y[i], average, y[i], col=adjustcolor(color, alpha.f=0.70), lwd=1.70)
            offset <- seq(0.10, 0.27, length.out=length(values))
            points(values, y[i]+offset, pch=16, col=adjustcolor(color, alpha.f=0.59), cex=0.70)
        } else {
            draw_distribution(values, y[i], color, horizontal=TRUE)
        }
        points(average, y[i], pch=if (stem) 21 else 23, cex=if (hero) 1.58 else 1.25,
               bg=color, col="white")
        text(labels_x, y[i], if (stem) sprintf("%.1f%%", 100*average) else sprintf("%.4f", average),
             adj=1, cex=0.80, col=if (hero) color else ink, font=if (hero) 2 else 1)
        mtext(labels[method], side=2, at=y[i], line=0.46, las=1, adj=1, cex=0.83, col=ink,
              font=if (hero) 2 else 1)
    }
    axis(1, at=ticks, labels=tick_labels, cex.axis=0.9)
    box(bty="l", col=ink)
    panel(letter, title)
}

save_figure("fig02_pattern_d_ablation", 11.35, 4.55, function() {
    layout(matrix(1:2, nrow=1), widths=c(1.07, 1))
    methods <- c("E0", "E4", "E5", "E6", "E7", "E8", "Shannon", "AL2CO variance")
    par(mar=c(3.6, 10.1, 2.8, 0.8))
    d_panel(methods, "medium_as_variable", c(-0.035,1.23), seq(0,1,0.25),
            paste0(seq(0,100,25), "%"), "Medium-rate sites classified as variable", "A", "Thresholded detection", stem=TRUE)
    par(mar=c(3.6, 10.1, 2.8, 1.0))
    d_panel(methods, "spearman", c(0.795,1.055), seq(0.8,1,0.05),
            sprintf("%.2f", seq(0.8,1,0.05)), "Spearman correlation", "B", "Three-rate ranking")
})

profile_panel <- function(seed, methods, limits, xticks, line_types, letter, title, legend_at) {
    values <- profiles[profiles$seed == seed, ]
    example <- values[values$method == methods[1], ]
    visible <- values[values$position >= limits[1] & values$position <= limits[2], ]
    ylimits <- range(visible$margin, 0)
    ylimits <- ylimits + diff(ylimits) * c(-0.09,0.21)
    yticks <- pretty(ylimits, 5)
    yticks <- yticks[yticks >= ylimits[1] & yticks <= ylimits[2]]
    plot(NA, xlim=limits, ylim=ylimits, xlab="Alignment position (1-based)",
         ylab="Score minus frozen threshold", axes=FALSE, xaxs="i", yaxs="i")
    runs <- rle(example$rate)
    ends <- cumsum(runs$lengths) + 0.5
    starts <- c(0.5, head(ends, -1))
    truth_colors <- c("0.06"="#E1EEE4", "0.32"="#F5EDD5", "1"="#FFFFFF")
    for (i in seq_along(ends)) {
        rect(starts[i], ylimits[1], ends[i], ylimits[2], col=truth_colors[as.character(runs$values[i])], border=NA)
    }
    axis_frame(xticks, yticks)
    abline(h=0, lty=3, col=muted, lwd=1.15)
    for (i in seq_along(methods)) {
        data <- values[values$method == methods[i], ]
        lines(data$position, data$margin, col=colors[methods[i]], lty=line_types[i], lwd=1.65)
    }
    legend(legend_at, methods, col=colors[methods], lwd=1.65, lty=line_types, bty="n",
           horiz=TRUE, cex=0.80, seg.len=2)
    panel(letter, title)
}

save_figure("fig03_position_profiles", 11.8, 4.04, function() {
    layout(matrix(1:2, nrow=1), widths=c(1,1.37))
    par(mar=c(3.6, 4.4, 2.8, 1.0))
    profile_panel(50124, c("E2", "E3", "E8"), c(105,136), c(105,110,115,120,125,130,135),
                  c(2,3,1), "A", "C | 3-residue motif", "topleft")
    profile_panel(50140, c("E5", "E7", "E8"), c(1,240), c(1,40,80,120,160,200,240),
                  c(2,3,1), "B", "D | 40% high-rate block", "topright")
})

save_figure("fig04_pattern_d_roc_pr", 10.8, 4.4, function() {
    layout(matrix(1:2, nrow=1))
    chosen <- c("E0", "E5", "E7", "E8")
    styles <- c(2,4,3,1)
    for (type in c("roc", "pr")) {
        par(mar=c(3.65, 4.3, 3.0, 0.75))
        is_roc <- type == "roc"
        xlim <- if (is_roc) c(0,0.10) else c(0,1)
        ylim <- c(0.50,1.01)
        xticks <- if (is_roc) seq(0,0.10,0.02) else seq(0,1,0.2)
        plot(NA, xlim=xlim, ylim=ylim, xlab=if (is_roc) "False-positive rate" else "Recall",
             ylab=if (is_roc) "True-positive rate" else "Precision", axes=FALSE, xaxs="i", yaxs="i")
        axis_frame(xticks, seq(0.5,1,0.1), if (is_roc) sprintf("%.2f", xticks) else xticks)
        if (!is_roc) {
            abline(h=curves$prevalence[1], col=muted, lty=2, lwd=1)
            text(0.03, curves$prevalence[1]+0.016, "Variable-site prevalence", adj=0,
                 cex=0.77, col=muted)
        }
        values <- numeric(length(chosen))
        for (i in seq_along(chosen)) {
            selected <- curves[curves$method == chosen[i], ]
            if (is_roc) {
                lines(selected$fpr, selected$recall, col=colors[chosen[i]], lty=styles[i], lwd=1.75)
                values[i] <- selected$auroc[1]
            } else {
                lines(selected$recall, selected$precision, col=colors[chosen[i]], lty=styles[i], lwd=1.75, type="S")
                values[i] <- selected$ap[1]
            }
        }
        legend("bottomright", paste(chosen, sprintf("%.4f", values)),
               title=if (is_roc) "Full AUROC" else "Average precision",
               col=colors[chosen], lty=styles, lwd=1.65, bty="n", cex=0.83, seg.len=2.2)
        panel(if (is_roc) "A" else "B", if (is_roc) "ROC at low false-positive rates" else "Precision-recall")
    }
})

save_figure("fig05_reference_gap_sweep", 9.2, 3.65, function() {
    layout(matrix(1:2, nrow=1))
    chosen <- c("E5", "E6", "E8")
    for (metric in c("scored_fraction", "boundary_f1")) {
        coverage <- metric == "scored_fraction"
        par(mar=c(3.6, 4.4, 2.9, 0.9))
        plot(NA, xlim=c(0.07,0.53), ylim=c(-0.04,1.05), xlab="Nominal gap fraction",
             ylab=if (coverage) "Scored fraction" else "Boundary F1", axes=FALSE, xaxs="i", yaxs="i")
        axis_frame(c(0.1,0.3,0.5), seq(0,1,0.25), c("10%","30%","50%"))
        for (i in seq_along(if (coverage) "E8" else chosen)) {
            method <- (if (coverage) "E8" else chosen)[i]
            selected <- gap[gap$method == method & gap$metric == metric, ]
            selected <- selected[order(selected$nominal_gap), ]
            x <- selected$nominal_gap + if (coverage) 0 else c(-0.013,0,0.013)[i]
            spread <- selected$sd > 0
            arrows(x[spread], pmax(0, selected$mean-selected$sd)[spread], x[spread],
                   pmin(1, selected$mean+selected$sd)[spread],
                   col=colors[method], lwd=1.3, length=0.045, code=3, angle=90)
            lines(x, selected$mean, type="b", lwd=1.9, col=colors[method],
                  pch=c(19,17,15)[i], lty=c(1,2,1)[i])
        }
        if (!coverage) legend("bottomleft", chosen, col=colors[chosen], lwd=1.8, bty="n",
                              pch=c(19,17,15), lty=c(1,2,1), horiz=TRUE, cex=0.85)
        panel(if (coverage) "A" else "B", if (coverage) "Usable columns" else "Boundary recovery")
    }
})
