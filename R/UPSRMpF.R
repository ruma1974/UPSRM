# _"/media/rune/joachime/RscripterCool/Rprojects/rmpackges/UPSRMp.R"



#' loads MaxQuant results
#' @param Dir (path to MaxQuant result)
#' ---
#' @export
loadMaxQuant <- function(Dir){
#
#
if (grepl("/",Dir)){
DFub=read.delim(paste0(Dir,"/GlyGly (K)Sites.txt"), sep="\t", stringsAsFactors=F,header=T, na.strings="")
DFpro=read.delim(paste0(Dir,"/proteinGroups.txt"), sep="\t", stringsAsFactors=F,header=T, na.strings="")
} else {
DFub=read.delim(paste0(Dir,"\\GlyGly (K)Sites.txt"), sep="\t", stringsAsFactors=F,header=T, na.strings="")
DFpro=read.delim(paste0(Dir,"\\proteinGroups.txt"), sep="\t", stringsAsFactors=F,header=T, na.strings="")
}
res=list()
res$DFub=DFub
res$DFpro=DFpro
return(res) # invisible(res)
}

#' return the names for possible UPS groups (to use with getQdata)
#' ---
#' @export
showUPSgroups <- function(){
return(names(UBL$UPSL)) # invisible(res)
}

#' Return gene names for the specific UPS group
#' @param name_s (name of UPS group - see showUPSgroups)
#' ---
#' @export
showGenes <- function(name_s){
return(UBL$UPSL[[name_s]]) # invisible(res)
}

#' return quantitative LFQ data for UPS group
#' @param res_l (loaded MaxQuant data)
#' @param UPSname_s (name of UPS group - see showUPSgroups)
#' @return
#' @keywords
#' @examples
#' ---
#' @export
getQdata <- function(res_l,UPSname_s){
#
#
idxL=unlist(lapply(res_l$DFpro$Gene.names,function(x){any(unlist(strsplit(x, split=";")) %in% UBL$UPSL[[UPSname_s]],na.rm=TRUE)}))
iL=grep("LFQ.",names(res_l$DFpro))
res=res_l$DFpro[idxL,iL]
res$GN=res_l$DFpro$Gene.names[idxL]
#
DF=res
gnL=DF$GN
DF=DF[,-which("GN"==names(DF))]
rownames(DF)=gnL
mat=DF
m = as.data.frame(mat)
m$keyvalueXXX = rownames(mat)
m = plyr::ddply(m, "keyvalueXXX", plyr::numcolwise(mean))
idx = -which("keyvalueXXX" == names(m))
res = m[, idx]
res = as.matrix(res)
rownames(res) = m[, -idx]
return(res) # invisible(res)
}

#' Return possible ubiquitin and ubiquitin like modifiers - to use with getUbQdata
#' ---
#' @export
showUBls <- function(){
#
#
res=c("Ub","SUMO1","SUMO2","SUMO3","SUMO4","NEDD8","UBD","ATG12","GABARAP","GABARAPL1","GABARAPL2","MAP1LC3B2","MAP1LC3B","ISG15","UFM1")
return(res) # invisible(res)
}


#' Return quantitative data for a UbL
#' @param res_l (MaxQuant result)
#' @param UPSname_s (see showUBls)
#' @param Fac (factor - if set to null only quantitative data is returned)
#' @param comp (comparisons to make e.g. list( c("MG", "DMSO"), c("PR", "DMSO")))
#' @param filename (filename for PDF output)
#' ---
#' @export
getUbQdata <- function(res_l,UPSname_s,Fac=NULL,comp=NULL,filename=NULL){
#
#
if (UPSname_s=="Ub"){DFcur=UBL$DFubRef} # Ub
if (UPSname_s=="SUMO1"){DFcur=UBL$DFsumo1Ref} # SUMO1
if (UPSname_s=="SUMO2"){DFcur=UBL$DFsumo2Ref} # SUMO2
if (UPSname_s=="SUMO3"){DFcur=UBL$DFsumo3Ref} # SUMO3
if (UPSname_s=="SUMO4"){DFcur=UBL$DFsumo4Ref} # SUMO4
if (UPSname_s=="NEDD8"){DFcur=UBL$DFneed8Ref} # NEDD8
if (UPSname_s=="UBD"){DFcur=UBL$DFubdRef} # UBD
if (UPSname_s=="ATG12"){DFcur=UBL$DFatg12Ref} # ATG12
if (UPSname_s=="GABARAP"){DFcur=UBL$DFGABARAPref} # GABARAP
if (UPSname_s=="GABARAPL1"){DFcur=UBL$DFGABARAPL1ref} # 
if (UPSname_s=="GABARAPL2"){DFcur=UBL$DFGABARAPL2ref} #
if (UPSname_s=="MAP1LC3B2"){DFcur=UBL$DFMAP1LC3B2ref} #
if (UPSname_s=="MAP1LC3B"){DFcur=UBL$DFMAP1LC3B2ref} #
if (UPSname_s=="MAP1LC3A"){DFcur=UBL$MAP1LC3Aref} #
if (UPSname_s=="ISG15"){DFcur=UBL$DFISG15ref} #
if (UPSname_s=="UFM1"){DFcur=UBL$DFUFM1ref} #
#print(DFcur)
idxL=which(grepl("Intensity.",names(res_l$DFub)) & !grepl("___",names(res_l$DFub)))
res_l$DFub$seq=gsub(r"{\s*\([^\)]+\)}","",as.character(res_l$DFub$GlyGly..K..Probabilities))
iL=unlist(lapply(res_l$DFub$seq,function(x){x %in% DFcur[,c("seq")]}))
idxL=append(idxL,ncol(res_l$DFub))
resO=res_l$DFub[iL,idxL]	

resO=resO %>% left_join( DFcur, by = "seq")
M=resO[,1:(ncol(resO)-2)]

#x11()
#boxplot(as.numeric(M[1,])~as.factor(Fac))
#library("ggpubr")
if (!is.null(Fac)){
pdf(filename,w=10,h=10)

for (i in it(DFcur)){
ggDF=data.frame(samples=Fac,y=log2(1+as.numeric(M[i,])))
p <- ggpubr::ggboxplot(ggDF, x = "samples", y = "y",
                color = "samples", palette =rainbow(length(levels(as.factor(Fac)))),
                add = "jitter", shape = "samples",main=paste("Site",DFcur[i,c("Position")]))

p= p +  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
p= p + ggplot2::theme(plot.title = ggplot2::element_text(size = 20, face = "bold"))
p = p + ggplot2::labs(x = "Samples", y = bquote(~Log[2] ~ italic(Intensity)))
if (!is.null(comp)){
p=p + ggpubr::stat_compare_means(comparisons = comp)+ ggpubr::stat_compare_means(label.y = max(ggDF$y))# Add global p-value # Add pairwise comparisons p-value	,method = "t.test"
} # 
print(p)
} # i

dev.off()
} # is.null


return(resO) # invisible(res)
}
