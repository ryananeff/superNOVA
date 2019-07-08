plotGOVolcano <- function(supergo_file, diffexp_file, gene_cutoff = 0.05, num_plots=10){
  diffexp = read.table(diffexp_file,header = T,sep="\t")
  rownames(diffexp) = make.names(diffexp[,1],unique = T) #gene names should be in the first column and should be unique
  diffexp[,"logP"] = -1*log(diffexp$P.Value)

  supergo_res = read.table(supergo_file,header=T,sep="\t")
  supergo_res[,"logP"] = sapply(-1*log(supergo_res$pVal),FUN=function(x){min(x,50)})
  supergo_res = supergo_res[order(supergo_res[,"logP"],decreasing = T),]

  for (ix in 1:num_plots){
    top_mod = supergo_res[ix,]
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

    gene_pval = sapply(unlist(strsplit(as.character(top_mod$gene_pVal),split = '; ')),
                  FUN=function(x){min(-1*log(as.numeric(unlist(strsplit(as.character(x),split=':'))[2])),12)},USE.NAMES = F)
    names(gene_pval) = sapply(unlist(strsplit(as.character(top_mod$gene_pVal),split = '; ')),
                              FUN=function(x){unlist(strsplit(as.character(x),split=':'))[1]},USE.NAMES = F)

    gene_pval = gene_pval[names(gene_changecor)]

    genes = names(gene_pval)
    diffexp_row = diffexp[names(gene_pval),"logFC"]
    sizes = diffexp[names(gene_pval),"logP"]

    data = data.frame(x=gene_changecor,y=gene_pval,c=diffexp_row,s=sizes,name=names(gene_pval))

    ggplot(data,aes(x=x,y=y,color=c,size=s)) +
          geom_point() +
          scale_color_gradient2(low="blue",high="red",mid="white") +
          geom_text_repel(data = subset(data,y>3)[1:30,],aes(label=name),color="brown") +
          theme(text=element_text(size=14)) +
          labs(title=name,x="Mean correlation change",y="Gene level -logP",size="DEG -logP",color="DEG logFC") +
          geom_hline(yintercept=3,linetype="dashed",color="red") +
          geom_vline(xintercept=0,linetype="dashed",color="red")

    invisible(readline(prompt="Press [enter] to continue"))
  }
}
