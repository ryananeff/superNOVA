plotGOVolcano <- function(supergo_file, diffexp_file, gene_cutoff = 0.05, pathway_cutoff = 0.05,
                          num_plots=10,num_genes=20,
                          prefix="plotGO_",plt_size=5){

  options(stringsAsFactors = F)
  options(warn = 0)
  library("ggplot2")
  library("ggrepel")

  diffexp = read.table(diffexp_file,header = T,sep="\t")
  rownames(diffexp) = make.names(diffexp[,1],unique = T) #gene names should be in the first column and should be unique
  diffexp[,"logP"] = -1*log(as.numeric(diffexp$P.Value),base=10)

  supergo_res = read.table(supergo_file,header=T,sep="\t")
  supergo_res[,"logP"] = sapply(-1*log(as.numeric(supergo_res$pVal),base=10),FUN=function(x){min(x,20)})
  supergo_res = supergo_res[order(supergo_res[,"logP"],decreasing = T),]

  for (ix in 1:num_plots){
    top_mod = supergo_res[ix,]
    if (top_mod$pVal > pathway_cutoff){ next }
    name = top_mod$Module

    gene_changecor_up = sapply(unlist(strsplit(as.character(top_mod$Top_GOC),split = '; ')),
                            FUN=function(x){as.numeric(unlist(strsplit(as.character(x),split=':'))[2])},USE.NAMES = F)
    names(gene_changecor_up) = sapply(unlist(strsplit(as.character(top_mod$Top_GOC),split = '; ')),
                                   FUN=function(x){unlist(strsplit(as.character(x),split=':'))[1]},USE.NAMES = F)

    gene_changecor_dn = sapply(unlist(strsplit(as.character(top_mod$Top_LOC),split = '; ')),
                               FUN=function(x){as.numeric(unlist(strsplit(as.character(x),split=':'))[2])},USE.NAMES = F)
    names(gene_changecor_dn) = sapply(unlist(strsplit(as.character(top_mod$Top_LOC),split = '; ')),
                                      FUN=function(x){unlist(strsplit(as.character(x),split=':'))[1]},USE.NAMES = F)

    gene_changecor = c(gene_changecor_up,gene_changecor_dn)
    gene_changecor = gene_changecor[unique(names(gene_changecor))]

    gene_pval = sapply(unlist(strsplit(as.character(top_mod$gene_pVal),split = '; ')),
                  FUN=function(x){min(-1*log(as.numeric(unlist(strsplit(as.character(x),split=':'))[2]),base=10),7)},USE.NAMES = F)
    names(gene_pval) = sapply(unlist(strsplit(as.character(top_mod$gene_pVal),split = '; ')),
                              FUN=function(x){unlist(strsplit(as.character(x),split=':'))[1]},USE.NAMES = F)

    gene_pval = gene_pval[names(gene_changecor)]

    genes = names(gene_pval)
    diffexp_row = diffexp[names(gene_pval),"logFC"]
    sizes = diffexp[names(gene_pval),"logP"]

    data = data.frame(x=gene_changecor,y=gene_pval,c=diffexp_row,s=sizes,name=genes)

    pdf(paste0(prefix,name,".pdf"),width=plt_size*1.2,height=plt_size)
    g = ggplot(data,aes(x=x,y=y,color=c,size=s)) +
          xlim(-0.5,0.5) + ylim(-0.3,10) +
          geom_point() +
          scale_color_gradient2(low="blue",high="red",mid="white") +
          geom_text_repel(data = subset(data[order(gene_pval,decreasing = T)[1:num_genes],],y>-1*log(gene_cutoff,base=10)),
                          aes(label=name),
                          color="brown",
                          box.padding=0.7,
                          point.padding=0.2,
                          segment.color="gray") +
          theme(text=element_text(size=14),plot.title=element_text(size=10)) +
          labs(title=name,x="Mean correlation change",y="Gene level -logP",size="DEG -logP",color="DEG logFC") +
          geom_hline(yintercept=-1*log(gene_cutoff,base=10),linetype="dashed",color="red") +
          geom_vline(xintercept=0,linetype="dashed",color="red")
    print(g)
    dev.off()
    message(paste("Plotting",name,"..."))
    #invisible(readline(prompt="Press [enter] to continue"))
  }
}
