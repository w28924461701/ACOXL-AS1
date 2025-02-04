

#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("ggExtra")

library(ggplot2)
library(ggpubr)
library(ggExtra)

setwd("D:\\panCancer\\22.CIBERSORTcor")             #设置工作目录
pFilter=0.001         #设置p值过滤条件

#读取表达文件
exp=read.table("singleGeneExp.txt", header=T,sep="\t",row.names=1,check.names=F)
gene=colnames(exp)[1]
#读取免疫细胞浸润文件
TME=read.table("CIBERSORT.result.txt", header=T,sep="\t",row.names=1,check.names=F)
#去除正常样品
group=sapply(strsplit(row.names(exp),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
exp=exp[group==0,]
#样品取交集
sameSample=intersect(row.names(TME),row.names(exp))
TME=TME[sameSample,]
exp=exp[sameSample,]

#相关性检验
outTab=data.frame()
#按肿瘤类型循环
for(i in levels(exp[,"CancerType"])){
    exp1=exp[(exp[,"CancerType"]==i),]
    TME1=TME[(TME[,"CancerType"]==i),]
    y=as.numeric(exp1[,1])
    outVector=data.frame(i,gene)
	#按免疫细胞循环
	for(j in colnames(TME1)[1:22]){
		x=as.numeric(TME1[,j])
		if(sd(x)>0.01){
			df1=as.data.frame(cbind(x,y))
			corT=cor.test(x,y,method="spearman")
			cor=corT$estimate
			pValue=corT$p.value
			outVector=cbind(outVector,pValue)
			p1=ggplot(df1, aes(x, y)) + 
				xlab(j)+ylab(gene)+
				ggtitle(paste0("Cancer: ",i))+theme(title=element_text(size=10))+
			    geom_point()+ geom_smooth(method="lm") + theme_bw()+
			    stat_cor(method = 'spearman', aes(x =x, y =y))
		    p2=ggMarginal(p1, type = "density", xparams = list(fill = "orange"),yparams = list(fill = "blue"))
			if(pValue<pFilter){
				pdf(file=paste0("estimateCor.",i,"_",j,".pdf"),width=5,height=5)
				print(p2)
				dev.off()
			}
		}
		else{
			outVector=cbind(outVector,pValue=1)
		}
	}
	outTab=rbind(outTab,outVector)
}
colNames=c("CancerType","Gene",colnames(TME)[1:22])
colnames(outTab)=colNames
write.table(outTab,file="CIBERSORTcor.result.txt",sep="\t",row.names=F,quote=F)


