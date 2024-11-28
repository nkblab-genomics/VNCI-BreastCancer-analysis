#copy number segmentation to be obtained through ASCAT
#structural variation to be obtained through DELLY2 followed by filtering out imprecise SVs
#create a text file "mysamples.txt" with tumour sample IDs.
#Dependency: ShatterSeek, BioCircos
#Change the BioCircos parameters as required

library(ShatterSeek)
library(BioCircos)

runss = function(samplex) {
cnv_data = read.table(paste0("data/example/shattersheek_input_cnv/", samplex, "_sgatterSheek.cn.tsv"), header = T, sep = "\t")
sv_data = read.table(paste0("data/example/shattersheek_input_sv/", samplex, "_sv_shatterSheek.tsv"), header = T, sep = "\t")

SV_data <- SVs(chrom1=as.character(sv_data$chrom1), 
               pos1=as.numeric(sv_data$pos1),
               chrom2=as.character(sv_data$chrom2), 
               pos2=as.numeric(sv_data$pos2),
               SVtype=as.character(sv_data$SVtype), 
               strand1=as.character(sv_data$strand1),
               strand2=as.character(sv_data$strand2))

CN_data <- CNVsegs(chrom=as.character(cnv_data$chrom),
                   start=cnv_data$start,
                   end=cnv_data$end,
                   total_cn=cnv_data$CN)

chromothripsis <- shatterseek(
  SV.sample=SV_data,
  seg.sample=CN_data,
  genome="hg19")

write.table(chromothripsis@chromSummary, file = paste0("xx_", samplex, "_shatterSheek_out.tsv"), sep = "\t", quote = F)
}

samples = readLines("mysamples.txt")
for(i in samples){
  runss(i)
}

# combine output of resistant and sensitive groups in "all.resistant.tsv", and "all.sensetive.tsv".
# the above files must contain headers: sample, chrom, start, end, oscl_cn2.

resi_data1 = read.table("all.resistant.tsv", sep = '\t', header = T)
mysamples = unique(resi_data1$sample)
for(i in 2:20){
  resi_data = resi_data1 %>% filter(sample == mysamples[i]) %>% filter(oscl_cn2 < 5)
  snvChr = resi_data$chrom
  snvStart = resi_data$start
  snvEnd = resi_data$end
  snvValues = rep(i,dim(resi_data)[1])
  tracks = tracks + BioCircosCNVTrack('cnv_track', as.character(snvChr), snvStart, snvEnd, snvValues, 
                                      color = "#fcae91", range = c(0,25))
}

sens_data1 = read.table("all.sensetive.tsv", sep = '\t', header = T)
mysamples = unique(sens_data1$sample)
for(i in 1:20){
  sens_data = sens_data1 %>% filter(sample == mysamples[i]) %>% filter(oscl_cn2 < 5)
  snvChr = sens_data$chrom
  snvStart = sens_data$start
  snvEnd = sens_data$end
  snvValues = rep(i,dim(resi_data)[1])
  tracks = tracks + BioCircosCNVTrack('cnv_track', as.character(snvChr), snvStart, snvEnd, snvValues+20, 
                                      color = "#74a9cf", range = c(0,25))
}

tracks = tracks + BioCircosBackgroundTrack("arcs_background", colors = "red", fillColors = "white",
                                           borderColors = "red",
                                           minRadius = 0.63, maxRadius = 0.93 )
BioCircos(tracks, genomeTicksDisplay = F, genomeLabelDy = 0, yChr = F, genomeFillColor = rep("lightgrey",23), genome = "hg19")
