#circlize test
library(stringr) 
library(circlize)      
library(ComplexHeatmap)     
library(grid)    
library(extrafont)
loadfonts()

#circos.genomicInitialize(seq_stat)

circos.clear()
#circle_size = unit(1, "snpc")
#circos.par(gap.degree = 3)
#circos.genomicInitialize(seq_stat, plotType = 'axis')



pdf("./TJ1_SVs.distribution_plot.pdf", width = 14, height = 8)

circos.initializeWithIdeogram(species = "hg38", chromosome.index = paste0("chr", c(1:22)))

SV_message_data = read.table("./TJ1_SVs.prepare_for_circlize.txt",header=T,stringsAsFactors = FALSE)
depth_base <- as.data.frame(SV_message_data)


#######绘制SV
#########DEL
DEL <- subset(depth_base,SV_type == "DEL")
circos.trackHist(DEL$chr,track.height = 0.13,lwd=2.5,
                x=DEL$start,col = "#F6511D", border = "#F6511D", 
                bg.border = NA,bin.size = 1000000)


#######INS
INS <- subset(depth_base,SV_type == "INS")
circos.trackHist(INS$chr,track.height = 0.13,lwd=2.5,
                x=INS$start,col = "#C48900", border = "#C48900", 
                bg.border = NA,bin.size = 1000000)

#######DUP
DUP <- subset(depth_base,SV_type == "DUP")
circos.trackHist(INS$chr,track.height = 0.13,lwd=2.5,
                x=INS$start,col = "#007DB2", border = "#007DB2", 
                bg.border = NA,bin.size = 1000000)

#######INV
INV <- subset(depth_base,SV_type == "INV")
circos.trackHist(INS$chr,track.height = 0.13,lwd=2.5,
                x=INS$start,col = "#659300", border = "#659300", 
                bg.border = NA,bin.size = 1000000)


sv_legend <- Legend(
    at = c(1, 2, 3, 4), labels = c(' DEL (7,932)', ' INS (12,561) ', ' DUP (321)',' INV (76)'),
    labels_gp = gpar(fontsize = 14), grid_height = unit(0.7, 'cm'), grid_width = unit(0.7, 'cm'), 
    type = c('points', 'points', 'lines', 'lines'), pch = NA, background = c(NA, NA, NA, NA),
    legend_gp = gpar(col = c("#F6511D", "#C48900", "#007DB2",  "#659300"), lwd = 5),gap = 0.3,
    title = 'SV type', title_gp = gpar(fontsize = 16))

y_coord <- 0.5
x_coord <- 0.51
pushViewport(viewport(x = x_coord, y = y_coord))
grid.draw(sv_legend)
upViewport()
dev.off()


