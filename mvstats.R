## Open Source written by Kevin McGarigal of University of Massachuestts Amerst


"box.plots" <-
function(x,var='',by='',save.plot=FALSE,
	col='blue',las=1,...){

oldpar<-par(no.readonly=TRUE)

if(!var==''){
	y<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	y<-as.data.frame(y)
	}
else{y<-as.data.frame(x)}
	
if(by==''){ #box-and-whisker w/o groups
	par(mfrow=c(1,1),mar=c(5,5,4,2)) #graphics settings
	for(i in 1:ncol(y)){ #loop thru variables
		boxplot(y[i],ylab=names(y[i]),col=col,las=las,
		main=paste('Box-and-Whisker Plot of',names(y[i]),sep=' '),...)
		if(save.plot==TRUE){
			dev.print(jpeg,file=paste('box.',names(y[i]),'.jpg',sep=''),width=600,height=800)
			} #end save
		if(!i==ncol(y)) {readline("Press return for next plot ")}
		} #end loop thru variables
	} #end bw w/o groups

else{ #box-and-whisker w/ groups
	n<-by.names(x,by) #create by variable
	y<-cbind(n,y) #bind with selected variables
	par(mfrow=c(1,1),mar=c(5,5,4,2)) #graphics settings
	for(i in 3:ncol(y)){ #loop thru variables
		boxplot(y[[i]]~y[[2]],col=col,las=las,
		ylab=names(y[i]),main=paste('Box-and-Whisker Plot of',names(y[i]),sep=' '),...)
		if(save.plot==TRUE){
			dev.print(jpeg,file=paste('box.',names(y[i]),'.jpg',sep=''),width=600,height=800)
			} #end save
		if(!i==ncol(y)) {readline("Press return for next plot ")}
		} #end loop thru variables
	} #end bw w/groups
par(oldpar)
} #end function

"by.names" <-
function(infile,by=names(infile)){

x <- unique(infile[,by, drop=FALSE])
z <- as.character(x[,1])
if(1 < length(by))
	for(i in 2:length(by)) {
		z <- paste(z,as.character(x[,i]), sep='.')
	}
z <- data.frame(z)
names(z) <- '..key'
t <- paste(by,collapse='.')
z <- cbind(z,x,..id = 1:dim(z)[1])
infile <- cbind(1:dim(infile)[1],infile); names(infile)[1] <- '..seq'
z <- merge(infile[,c(by,'..seq')],z,by=by)
z <- z[sort(z$..seq,index.return=TRUE)$ix,]
row.names(z) <- z$..seq
z <- z[,c('..id','..key')]
names(z) <- c('id',t)
return(z)
}

"class.monte" <-
function(y,grouping='',prop=.5,type='qda',perm=1000,prior='',...){

y<-as.data.frame(y)
if(!grouping==''){
	grp<-as.factor(grouping)
	}
else{
	grp<-as.factor(y[,1]) #groups assumed to be in first column
	y<-y[,-1]
	}

G<-length(levels(grp))
z<-matrix(0,perm,G+2) #create blank matrix

for(i in 1:perm){
	ran.split(y,grp,prop)
	if(type=='qda'){
		y.da<-qda(calibrate,grouping=grp.cal,...)
		}
	else if(type=='lda'){
		y.da<-lda(calibrate,grouping=grp.cal,...)
		}
	y.pred<-predict(y.da,newdata=validate)
	y.table<-table(grp.val,y.pred$class)
	for(j in 1:G){
		z[i,j]<-y.table[j,j]/sum(y.table[j,])
		}
		z[i,G+1]<-sum(diag(y.table))/sum(y.table)
		if(!prior==''){
			z[i,G+2]<-tau(y.table,prior=prior)
			}		
		else{
			z[i,G+2]<-cohen.kappa(y.table)
			}
	}

z2<-as.data.frame(apply(z,2,function(x){ #calculate stats
	z2<-c(quantile(x,probs=c(0,.05,.5,.95,1),na.rm=TRUE),
    mean(x,na.rm=TRUE))
    names(z2)<-c('minimum','5th percentile','median','95th percentile','maximum','mean')
	z2<-round(z2,3) #round elements to 3 decimal places
	}))

if(!prior==''){
	colnames(z2)<-c(levels(grp),'Total','Tau')
	}
else{
	colnames(z2)<-c(levels(grp),'Total','Kappa')
	}

cat('Split-sample Cross-Validation Classification Summary\n')
cat('Correct Classification Rate:\n')

return(z2)	
}

"clus.composite" <-
function(x,grp){

groups<-as.factor(grp)
groups<-levels(groups)

#compute mean by cluster
x.by<-by(x,grp,mean)
	
#create composite clusters
z<-matrix(0,length(x),length(groups)) #create blank matrix
for(i in 1:length(groups)){
	z[,i]<-x.by[[i]]
	}
z<-as.data.frame(z)
rownames(z)<-colnames(x)
colnames(z)<-groups
z<-t(z)
return(z)
}

"clus.stats" <-
function(x,grp){

cv<<-function(x,na.rm) sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)*100 
x<-as.data.frame(x)
groups<-as.factor(grp)
groups<-levels(groups)

#compute kruskal-wallis rank sum test p-value
p.value<-rep(0,length(x))
for(i in 1:length(x)){
	p<-kruskal.test(x[,i],grp)
	p.value[i]<-p$p.value
	}
p.value<-round(p.value,3)

#compute number of obs. per cluster
cl.nobs<-rep(0,length(groups))
for(i in 1:length(groups)){
	cl.nobs[i]<-sum(grp==i)
	}

#compute stats by cluster
x.mean<-by(x,grp,mean)
x.cv<-by(x,grp,cv)

#create summary table
cl.mean<-matrix(0,length(x),length(groups)) #create blank matrix
cl.cv<-matrix(0,length(x),length(groups)) #create blank matrix
for(i in 1:length(groups)){
	cl.mean[,i]<-round(x.mean[[i]],3)
	cl.cv[,i]<-round(x.cv[[i]],0)
	}
cl.mean<-as.data.frame(cl.mean)
rownames(cl.mean)<-colnames(x)
colnames(cl.mean)<-groups
cl.mean<-cbind(cl.mean,p.value)

cl.cv<-as.data.frame(cl.cv)
rownames(cl.cv)<-colnames(x)
colnames(cl.cv)<-groups

z<-list(cl.nobs,cl.mean,cl.cv)
names(z)<-c('Cluster.nobs','Cluster.mean','Cluster.cv')
return(z)
}

"cohen.kappa" <-
function(y){

N<-sum(y)
ccr<-sum(diag(y))/sum(y)

p<-apply(y,1,sum)/N
q<-apply(y,2,sum)/N
num<-ccr-sum(p*q)
den<-1-sum(p*q)
kappa<-num/den
kappa[kappa<0]<-0
return(kappa)		
}

"contrast.matrix" <-
function(grp){

grp<-as.factor(grp)
N<-length(grp)
grp.mat<-matrix(1,nrow=N,ncol=N)
matched<-function(irow,icol,grp){
	    grp[irow]==grp[icol]
		}
irow<-row(matrix(nrow=N,ncol=N))
icol<-col(matrix(nrow=N,ncol=N))
grp.mat[matched(irow,icol,grp)]<-0
z<-as.dist(grp.mat)
return(z)
}

"cov.test" <-
function(x,groups,var='',method='bartlett',...){

if(!var==''){
	y<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	y<-as.data.frame(y)
	}
else{y<-as.data.frame(x)}

#create summary table
z<-matrix(0,ncol(y),2) #create blank matrix
for(i in 1:ncol(y)){ #loop thru variables
	if(method=='bartlett'){
		temp<-bartlett.test(y[,i],groups,...)
		}
	else if(method=='fligner'){
		temp<-fligner.test(y[,i],groups,...)
		}	
	z[i,1]<-temp$statistic
	z[i,2]<-temp$p.value
	}
z<-round(z,3)
rownames(z)<-colnames(y)

if(method=='bartlett'){
	colnames(z)<-c('Bartletts K-squared','p-value')
	cat('Bartlett Test of Homogeneity of Variances:\n')
	}
else if(method=='fligner'){
	colnames(z)<-c('Median chi-squared','p-value')
	cat('Fligner-Killeen Test of Homogeneity of Variances:\n')
	}

return(z)
}

"data.dist" <-
function(x,method,var='',cor.method='pearson',
	outfile='',binary=FALSE,diag=FALSE,upper=FALSE,na.rm=TRUE,...){

library(vegan) #load vegan library
library(MASS) #load MASS library

if(!var==''){
	y<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	}
else{y<-x}

#ecological distance calculations
if(method=='correlation'){ #compute correlation distance
	y<-t(y) #transpose dataset
	z<-as.dist((1-cor(y,method=cor.method,use='complete.obs',...))/2) #compute distance
	z<-round(z,3) #round results to 3 decimal places
	if(!outfile==''){ #save outfile
		zz<-as.matrix(z)
		write.table(zz,file=paste(outfile,'.csv',sep=''),quote=FALSE,sep=',')
		} #end save outfile
	} #end correlation distance
else { #compute all other distance methods
	z<-vegdist(y,method=method,binary=binary,diag=diag,
		upper=upper,na.rm=na.rm,...) #vegdist function
	z<-round(z,3) #round results to 3 decimal places
	if(!outfile==''){ #save outfile
		zz<-as.matrix(z)
		write.table(zz,file=paste(outfile,'.csv',sep=''),quote=FALSE,sep=',')
		} #end save outfile
	} #end all other distance methods
return(z)
} #end function

"data.stand" <-
function(x,method,var='',margin='column',
	outfile='',plot=TRUE,save.plot=FALSE,na.rm=TRUE,
	col.hist='blue',col.line='black',las=1,lab=c(5,5,4),...){

#things to do: add option for standardizing within groups

library(vegan) #load vegan library

if(plot==TRUE){
	old.par<-par(no.readonly=TRUE)
	}

x<-as.data.frame(x)

if(!var==''){
	y1<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	y2<-subset(x,select=-eval(parse(text=var))) #select remaining variables
	t1<-y1 #copy to work file for transformations
	}
else{
	y1<-x #original variables
	t1<-x #copy to work file for transformations
	}

#paired histogram function for comparing raw and standardized variable
paired.hist<-function(raw=y1,trans=t1){
	for(i in 1:ncol(raw)){ #loop thru selected variables
		par(mfrow=c(2,1),mar=c(5,5,4,2)) #graphics settings
		hist(raw[[i]],prob=TRUE,
		col=col.hist,las=las,lab=lab,
		xaxs='i',yaxs='i',xlab=names(raw[i]),
		main=paste('Histogram of',names(raw[i]),sep=' '),...)
		par(new=TRUE)
		plot(density(raw[[i]],na.rm=TRUE),
		col=col.line,las=las,lab=lab,
		xaxs='i',yaxs='i',axes=FALSE,xlab='',ylab='',main='',...)
		if(method=='wisconsin'){ #wisconsin plot
			hist(trans[[i]],prob=TRUE,
			col=col.hist,las=las,lab=lab,
			xaxs='i',yaxs='i',xlab=names(trans[i]),
			main=paste(method,'standardization',sep=' '),...)
			par(new=TRUE)
			plot(density(trans[[i]],na.rm=TRUE),
			col=col.line,las=las,lab=lab,
			xaxs='i',yaxs='i',axes=FALSE,xlab='',ylab='',main='',...)
			} #end wisconsin plot
		else{ #all other standardizations
			hist(trans[[i]],prob=TRUE,
			col=col.hist,las=las,lab=lab,
			xaxs='i',yaxs='i',xlab=names(trans[i]),
			main=paste(margin,method,'standardization',sep=' '),...)
			par(new=TRUE)
			plot(density(trans[[i]],na.rm=TRUE),
			col=col.line,las=las,lab=lab,
			xaxs='i',yaxs='i',axes=FALSE,xlab='',ylab='',main='',...)
			} #end all other plots
		if(save.plot==TRUE){ #save plot to file
			dev.print(jpeg,file=paste('shist.',names(raw[i]),'.jpg',sep=''),width=800,height=600)
			} #end save	plot
		readline("Press return for next plot ")
		} #end loop thru variables
	} #end paired histogram function

#standardizations
if (method=='wisconsin'){ #wisconsin standardization
	t1<-decostand(t1,method='max',na.rm=TRUE,...) #column max standardization
	t1<-decostand(t1,method='total',na.rm=TRUE,...) #row total standardization
	t1<-round(t1,3) #round results to 3 decimal places
	if(plot==TRUE){ #plot paired histograms
		paired.hist(y1,t1)
		} #end plot paired histograms
	if(!var==''){
		z<-cbind(y2,t1) #bind result to remaining variables
		}
	else{z<-t1}
	if(!outfile==''){ #save outfile
		write.table(z,file=paste(outfile,'.csv',sep=''),quote=FALSE,row.names=FALSE,sep=',')
		} #end save outfile
	} #end wisconsin standardization
else if (method=='chi.square'){ #chi-square standardization
	t1<-decostand(t1,method=method,MARGIN=1,na.rm=na.rm,...) #decostand function
	t1<-round(t1,3) #round results to 3 decimal places
	if(plot==TRUE){ #plot paired histograms
		paired.hist(y1,t1)
		} #end plot paired histograms
	if(!var==''){
		z<-cbind(y2,t1) #bind result to remaining variables
		}
	else{z<-t1}
	if(!outfile==''){ #save outfile
		write.table(z,file=paste(outfile,'.csv',sep=''),quote=FALSE,row.names=FALSE,sep=',')
		} #end save outfile
	} #end chi.square standardization
else { #all other standardizations
	if (margin=='row') {MARGIN=1} #assign MARGIN parameter for decostand(vegan)
	else if (margin=='column') {MARGIN=2} #assign MARGIN parameter for decostand(vegan)
	t1<-decostand(t1,method=method,MARGIN=MARGIN,na.rm=na.rm,...) #decostand function
	t1<-round(t1,3) #round results to 3 decimal places
	if(plot==TRUE){ #plot paired histograms
		paired.hist(y1,t1)
		} #end plot paired histograms
	if(!var==''){
		z<-cbind(y2,t1) #bind result to remaining variables
		}
	else{z<-t1}
	if(!outfile==''){ #save outfile
		write.table(z,file=paste(outfile,'.csv',sep=''),quote=FALSE,row.names=FALSE,sep=',')
		} #end save outfile
	} #end all other standardizations

if(plot==TRUE){
	par(old.par)
	}
return(z)
} #end function

"data.trans" <-
function(x,method,var='',exp=1,outfile='',
	plot=TRUE,save.plot=FALSE,col.hist='blue',col.line='black',
	las=1,lab=c(5,5,4),...){

if(plot==TRUE){
	old.par<-par(no.readonly=TRUE)
	}

if(!var==''){
	y1<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	y2<-subset(x,select=-eval(parse(text=var))) #select remaining variables
	t1<-y1 #copy to work file for transformations
	}
else{
	y1<-x #original variables
	t1<-x #copy to work file for transformations
	}

#paired histogram function for comparing raw and transformed variable
paired.hist<-function(raw,trans){
	for(i in 1:ncol(raw)){ #loop thru selected variables
		par(mfrow=c(2,1),mar=c(5,5,4,2)) #graphics settings
		hist(raw[[i]],prob=TRUE,col=col.hist,las=las,lab=lab,
			xaxs='i',yaxs='i',xlab=names(raw[i]),
			main=paste('Histogram of',names(raw[i]),sep=' '),...)
		par(new=TRUE)
		plot(density(raw[[i]],na.rm=TRUE),
			xaxs='i',yaxs='i',col=col.line,las=las,lab=lab,
			axes=FALSE,xlab='',ylab='',main='',...)
		if(method=='power'){ #plot title for power transform
			hist(trans[[i]],prob=TRUE,col=col.hist,las=las,lab=lab,
				xaxs='i',yaxs='i',xlab=names(trans[i]),
				main=paste(method,'(x^',exp,')','transformation',sep=' '),...)
			} #end plot title for power transform
		else{ #plot title for other transforms
			hist(trans[[i]],prob=TRUE,col=col.hist,las=las,lab=lab,
				xaxs='i',yaxs='i',xlab=names(trans[i]),
				main=paste(method,'transformation',sep=' '),...)
			} #end plot title for other transforms
		par(new=TRUE)
		plot(density(trans[[i]],na.rm=TRUE),
			xaxs='i',yaxs='i',col=col.line,las=las,lab=lab,
			axes=FALSE,xlab='',ylab='',main='',...)
		if(save.plot==TRUE){
			dev.print(jpeg,file=paste('thist.',method,'.',names(raw[i]),'.jpg',sep=''),width=800,height=600)
			} #end save	
		readline("Press return for next plot ")
		} #end loop thru variables
	} #end paired histogram function

#transformations
if(method=='power'){ #power transformations
	if(exp==0){ #binary transformation
		t1[t1!=0]<-1 #compute binary
		t1<-as.data.frame(round(t1,3)) #round to 3 decimal places
		if(plot==TRUE){
			paired.hist(y1,t1)
			}
		if(!var==''){
			z<-cbind(y2,t1) #bind transformed variables with remaining
			}
		else{z<-t1}
		if(!outfile==''){ #save outfile
			write.csv(z,file=outfile,quote=FALSE,row.names=FALSE)
			} #end save outfile
		} #end binary transformation
	else{ #other power transformations
		t1<-t1^exp #compute power transformation
		t1<-as.data.frame(round(t1,3)) #round to 3 decimal places
		if(plot==TRUE){
			paired.hist(y1,t1)
			}
		if(!var==''){
			z<-cbind(y2,t1) #bind transformed variables with remaining
			}
		else{z<-t1}
		if(!outfile==''){ #save outfile
			write.table(z,file=paste(outfile,'.csv',sep=''),quote=FALSE,row.names=FALSE,sep=',')
			} #end save outfile
		} #end other power transformations
	} #end power transformation

else if(method=='log'){ #log10 transformation
	q<-min(t1[t1!=0],na.rm=TRUE) #find minimum non-zero value
	c<-trunc(log10(q)) #order of magnitude constant
	d<-10^(c) #decimal constant (inverse log10 of c
	t1<-log10(t1+d)-c #compute log10 transformation
	t1<-as.data.frame(round(t1,3)) #round to 3 decimal places
	if(plot==TRUE){
		paired.hist(y1,t1)
		}
	if(!var==''){
		z<-cbind(y2,t1) #bind transformed variables with remaining
		}
	else{z<-t1}
	if(!outfile==''){ #save outfile
		write.table(z,file=paste(outfile,'.csv',sep=''),quote=FALSE,row.names=FALSE,sep=',')
		} #end save outfile
	} #end log transformation
	
else{ #arcsine squareroot transformation
	t1<-(2/pi)*asin(sqrt(t1)) #compute arcsin square root
	t1<-as.data.frame(round(t1,3)) #round to 3 decimal places
	if(plot==TRUE){
		paired.hist(y1,t1)
		}
	if(!var==''){
		z<-cbind(y2,t1) #bind transformed variables with remaining
		}
	else{z<-t1}
	if(!outfile==''){ #save outfile
		write.table(z,file=paste(outfile,'.csv',sep=''),quote=FALSE,row.names=FALSE,sep=',')
		} #end save outfile
	} #end arcsine squareroot transformation
	
if(plot==TRUE){
	par(old.par)
	}

return(z)
} #end function	

"dist.plots" <-
function(x,groups,distance='euclidean',na.rm=TRUE,
	col='blue',col.line='red',las=1,...){

#create vector of dissimilarities for plotting
if(inherits(x,'dist')) 
    dmat<-x
else if(is.matrix(x)&&nrow(x)==ncol(x)&&all(x[lower.tri(x)]==t(x)[lower.tri(x)])){
    dmat<-x
    attr(dmat,'method')<-'user supplied square matrix'
    }
else dmat<-vegdist(x,method=distance)
distance<-attr(dmat,'method')
dmat<-as.matrix(dmat)
diag(dmat)<-NA
x<-as.dist(dmat)
groups<-as.factor(groups)
matched<-function(irow,icol,groups){
    groups[irow]==groups[icol]
	}
diss.vec<-as.vector(x)
N<-attributes(x)$Size
irow<-as.vector(as.dist(row(matrix(nrow=N,ncol=N))))
icol<-as.vector(as.dist(col(matrix(nrow=N,ncol=N))))
within<-matched(irow,icol,groups)
class.vec<-rep("Between",length(diss.vec))
take<-as.numeric(irow[within])
class.vec[within]<-levels(groups)[groups[take]]
class.vec<-factor(class.vec,levels=c("Between",levels(groups)))

#boxplot of between and within-group dissimilarities
boxplot(diss.vec~class.vec,notch=TRUE,varwidth=TRUE,col=col,
	ylab='Dissimilarity',main='Between- and Within-group Dissimilarities',...)
readline("Press return for next plot ")

#histograms of w/i group dissimilarities
grp<-levels(groups) #create vector of group names
s<-floor(length(grp)^.5) #create multi-figure dimensions
s<-c(s,ceiling(length(grp)/s)) #create multi-figure dimensions
old.par<-par(no.readonly=TRUE)
par(mfrow=s,mar=c(5,5,4,2)) #graphics settings
for(j in 1:length(grp)){ #loop thru groups
	z<-diss.vec[class.vec==grp[j]] #select records for group j
	hist(z,prob=TRUE,col=col,las=las,xlim=range(diss.vec),
		xaxs='i',yaxs='i',xlab=grp[j],main='',...)
	par(new=TRUE)
	plot(density(z,na.rm=na.rm),xaxs='i',yaxs='i',
		col=col.line,las=las,xlim=range(diss.vec),
		axes=FALSE,xlab='',ylab='',main='',...)
	} #end loop thru groups
par(old.par)
readline("Press return for next plot ")

#qqnorm plots of w/i group dissimilarities
s<-floor(length(grp)^.5) #create multi-figure dimensions
s<-c(s,ceiling(length(grp)/s)) #create multi-figure dimensions
old.par<-par(no.readonly=TRUE)
par(mfrow=s,mar=c(5,5,4,2)) #graphics settings
for(j in 1:length(grp)){ #loop thru groups
	z<-diss.vec[class.vec==grp[j]] #select records for group j
	qqnorm(z,datax=TRUE,col=col,las=las,main=grp[j],
		ylab='Sample quantiles',...)
	z.IQR<-IQR(z,na.rm=na.rm)
	if(z.IQR>0){ #add qqline
		qqline(z,datax=TRUE,col=col.line,...)
		} #end qqline
	} #end loop thru groups
par(old.par)

}

"drop.var" <-
function(x,var='',outfile='',cv=0,minpo=0,minfo=0,
	maxpo=100,maxfo=nrow(x),pct.missing=100){

if(!var==''){
	y1<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	y2<-subset(x,select=-eval(parse(text=var))) #select remaining variables
	}
else{
	y1<-x #original variables
	}

#statistical functions
cv<<-function(x,na.rm) sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)*100 
freq.occur<<-function(x,na.rm) sum(!x==0,na.rm=TRUE)
pct.occur<<-function(x,na.rm) sum(!x==0,na.rm=TRUE)/length(x)*100
pct.missing<<-function(x,na.rm) sum(is.na(x))/length(x)*100 

#delete offending variables
z<-as.matrix(apply(y1,2,cv)<cv | 
apply(y1,2,pct.occur)>maxpo |
apply(y1,2,pct.occur)<minpo |
apply(y1,2,freq.occur)>maxfo |
apply(y1,2,freq.occur)<minfo | 
apply(y1,2,pct.missing)>pct.missing)
z<-y1[,z[,1]=="FALSE",,drop=FALSE] 
if(!var==''){
	z<-cbind(y2,z) #bind reduced set w/ deselected variables
	}
if(!outfile==''){ #save outfile
	write.table(z,file=paste(outfile,'.csv',sep=''),quote=FALSE,row.names=FALSE,sep=',')
	} #end save outfile
return(z)
} #end function

"ecdf.plots" <-
function(x,var='',by='',save.plot=FALSE,
	col.point='blue',col.line='red',las=1,lab=c(10,10,3),...){

old.par<-par(no.readonly=TRUE)

if(!var==''){
	y<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	y<-as.data.frame(y)
	}
else{y<-as.data.frame(x)}

if(by==''){ #ecdf w/o groups
	for(i in 1:ncol(y)){ #loop thru variables
		plot(sort(y[[i]]),type='o',col=col.point,lab=lab,las=las,
		yaxs='i',xaxs='i',
		xlab='Cumulative Number of Plots',ylab=names(y[i]),
		main='Empirical Cumulative Distribution Function',...)
		line(c(1,length(y[[i]])),range(y[[i]]))
		if(save.plot==TRUE){ #save plot
			dev.print(jpeg,file=paste('ecdf.',names(y[i]),'.jpg',sep=''),width=800,height=600)
			} #end save plot
		readline("Press return for next plot ")
		} #end loop thru variables
	} #end qqnorm w/o groups
	
else{ #ecdf plots w/ groups
	n<-by.names(x,by) #create by variable
	y<-cbind(n,y) #bind with selected variables
	groups<-levels(y[,2]) #create vector of group names
	s<-floor(length(groups)^.5) #create multi-figure dimension
	s<-c(s,ceiling(length(groups)/s)) #create multi-figure dimensions
	for(i in 3:ncol(y)){ #loop thru variables
		par(mfrow=s,mar=c(5,5,4,2))	#graphics settings	
		for(j in 1:length(groups)){ #loop thru groups
			z<-y[y[,2]==groups[j],] #select records for group j
			plot(sort(z[[i]]),type='o',col=col.point,lab=lab,las=las,
			yaxs='i',xaxs='i',
			xlab='Cum No of Plots',ylab=names(z[i]),
			main=groups[j],...)
			line(c(1,length(z[[i]])),range(z[[i]]))
			if(j==length(groups) & save.plot==TRUE){
				dev.print(jpeg,file=paste('ecdf.',names(z[i]),'.jpg',sep=''),width=800,height=600)
				} #end save
			} #end loop thru groups
			readline("Press return for next plot ")
		} #end loop thru variables
	} #end qqnorm w/ groups
par(old.par)
} #end function

"foa.plots" <-
function(x,var='',margin='column',outfile='',na.rm=TRUE,col.hist='blue',
	col.point='blue',col.line='red',col.den='black',las=1,lab=c(10,10,3),...){

old.par<-par(no.readonly=TRUE)

if(!var==''){
	y<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	}
else{y<-x}

z1<-apply(y>0,2,sum,na.rm=na.rm) #compute species frequency of occurrence
if(min(z1)==0) stop("all species must have non-zero sum")
z2<-100*z1/nrow(y) #compute species percent occurrence
z3<-log(z1) #compute log of species frequency of occurrences
z4<-apply(y,2,sum,na.rm=na.rm)/z1 #compute species mean abundance where present
r1<-apply(y>0,1,sum,na.rm=na.rm) #compute number of species per plot
r2<-apply(y,1,sum,na.rm=na.rm) #compute total abundance per plot

#summary table
if(margin=='column'){
	z<-as.data.frame(cbind(z1,z2,z3,z4))
	z<-round(z,3)
	names(z)<-c('spc.pres','spc.perc','spc.log','spc.mean')
	}
else{z<-as.data.frame(cbind(r1,r2))
	z<-round(z,3)
	names(z)<-c('plt.pres','plt.sum')
	}
if(!outfile==''){ #save outfile
	write.table(z,file=paste(outfile,'.csv',sep=''),quote=FALSE,sep=',')
	} #end save outfile
	
#plot1: cumulative distribution of species occurrence - raw scale
plot(sort(z1),type='o',col=col.point,lab=lab,las=las,
ylim=c(0,max(z1)),xlim=c(0,length(z1)),yaxs='i',xaxs='i',
xlab='Cumulative Number of Species',ylab='Frequency of Occurrence',
main='Cumulative Distribution of Species Occurrence',...)
abline(h=quantile(z1,.5),lty=1,col=col.line,...)
abline(h=quantile(z1,.05),lty=2,col=col.line,...)
abline(h=quantile(z1,.95),lty=2,col=col.line,...)
readline("Press return for next plot ")

#plot2: cumulative distribution of species occurrence - relativized
plot(sort(z2),type='o',col=col.point,lab=lab,las=las,
ylim=c(0,100),xlim=c(0,length(z2)),yaxs='i',xaxs='i',
xlab='Cumulative Number of Species',ylab='Percent Occurrence',
main='Cumulative Distribution of Species Relative Occurrence',...)
abline(h=50,lty=1,col=col.line,...)
abline(h=5,lty=2,col=col.line,...)
abline(h=95,lty=2,col=col.line,...)
readline("Press return for next plot ")

#plot3: histogram of species occurrence - raw scale
hist(z1,prob=TRUE,col=col.hist,las=las,lab=lab,
xaxs='i',yaxs='i',xlab='Species Occurrence',ylab='Frequency',
main='Histogram of Species Occurrence',...)
par(new=TRUE)
plot(density(z1,na.rm=na.rm),xaxs='i',yaxs='i',
col=col.den,las=las,lab=lab,
axes=FALSE,xlab='',ylab='',main='',...)
readline("Press return for next plot ")

#plot4: histogram of species occurrence - log transformation
hist(z3,prob=TRUE,col=col.hist,las=las,lab=lab,
xaxs='i',yaxs='i',xlab='Log(Species Occurrence)',ylab='Frequency',
main='Histogram of Log(Species Occurrence)',...)
par(new=TRUE)
plot(density(z3,na.rm=na.rm),xaxs='i',yaxs='i',
col=col.den,las=las,lab=lab,
axes=FALSE,xlab='',ylab='',main='',...)
readline("Press return for next plot ")

#plot5: cumulative distribution of species mean abundance - raw scale
plot(sort(z4),type='o',col=col.point,lab=lab,las=las,
ylim=c(0,max(z4)),xlim=c(0,length(z4)),yaxs='i',xaxs='i',
xlab='Cumulative Number of Species',ylab='Mean Abundance',
main='Cumulative Distribution of Species Mean Abundance',...)
abline(h=quantile(z4,.5),lty=1,col=col.line,...)
abline(h=quantile(z4,.05),lty=2,col=col.line,...)
abline(h=quantile(z4,.95),lty=2,col=col.line,...)
readline("Press return for next plot ")

#plot6: plot of frequency of occurrence versus mean abundance
plot(z1,z4,type='p',col=col.point,lab=lab,las=las,
yaxs='i',xaxs='i',
xlab='Frequency of Occurrence',ylab='Mean Abundance',
main='Species Occurrence vs Mean Abundance',...)
yorn<-readline("Do you want to identify individual species? Y/N :")
    if(yorn=='Y' || yorn=='y') 
        identify(z1,z4,names(y))
readline("Press return for next plot ")

#plot7: plot of frequency of occurrence versus log of mean abundance
plot(z1,z4,type='p',col=col.point,lab=lab,las=las,
yaxs='i',xaxs='i',log='y',
xlab='Frequency of Occurrence',ylab='Log(Mean Abundance)',
main='Species Occurrence vs Log(Mean Abundance)',...)
yorn<-readline("Do you want to identify individual species? Y/N :")
    if(yorn=='Y' || yorn=='y') 
        identify(z1,z4,names(y))
readline("Press return for next plot ")

#plot8: cumulative distribution of species per plot
plot(sort(r1),type='o',col=col.point,lab=lab,las=las,
ylim=c(0,max(r1)),xlim=c(0,length(r1)),yaxs='i',xaxs='i',
xlab='Cumulative Number of Plots',ylab='Plot Richness',
main='Cumulative Distribution of Plot Richness',...)
abline(h=quantile(r1,.5),lty=1,col=col.line,...)
abline(h=quantile(r1,.05),lty=2,col=col.line,...)
abline(h=quantile(r1,.95),lty=2,col=col.line,...)
readline("Press return for next plot ")

#plot9: cumulative distribution of total plot abundance
plot(sort(r2),type='o',col=col.point,lab=lab,las=las,
ylim=c(0,max(r2)),xlim=c(0,length(r2)),yaxs='i',xaxs='i',
xlab='Cumulative Number of Plots',ylab='Total Abundance',
main='Cumulative Distribution of Plot Total Abundance',...)
abline(h=quantile(r2,.5),lty=1,col=col.line,...)
abline(h=quantile(r2,.05),lty=2,col=col.line,...)
abline(h=quantile(r2,.95),lty=2,col=col.line,...)
readline("Press return for next plot ")

#plot10: plot of plot richness versus plot total abundance
plot(r1,r2,type='p',col=col.point,lab=lab,las=las,
yaxs='i',xaxs='i',
xlab='Plot Richness',ylab='Plot Total Abundance',
main='Plot Richness vs Total Abundance',...)
yorn<-readline("Do you want to identify individual plots? Y/N :")
    if(yorn=='Y' || yorn=='y') 
        identify(r1,r2,rownames(y))

par(old.par)
return(z)
} #end function

"hclus.cophenetic" <-
function(d,hclus,fit='lm',...){

old.par<-par(no.readonly=TRUE)
d.coph<-cophenetic(hclus)
r<-round(cor(d,d.coph),2)
plot(d,d.coph,xlab='Observed dissimilarity',
	ylab='Cophenetic dissimilarity',
	main=paste('Cophenetic Correlation ',
	'(',hclus$dist.method,', ',hclus$method,')',sep=''),...)
	text(max(d),min(d.coph),paste('Cophenetic correlation = ',r,sep=''),col='red',pos=2)
#	title(sub=paste('Cophenetic correlation = ',r,sep=''),col.sub='red',adj=0)
if(fit=='lm'){
	abline(lm(d.coph~d),col='blue')
	}
else if(fit=='rlm'){
	abline(rlm(d.coph~d),col='blue')
	}
else if(fit=='qls'){
	abline(lqs(d.coph~d),col='blue')
	}
	
par(old.par)
return(r)
}

"hclus.scree" <-
function(x,...){

old.par<-par(no.readonly=TRUE)
z1<-seq(length(x$labels)-1,1)
z<-as.data.frame(cbind(z1,x$height))
plot(z[,1],z[,2],type='o',lwd=1.5,pch=19,col='blue',
	ylab='Dissimilarity',xlab='Number of Clusters',
	main=paste('Scree Plot of Hierarchical Clustering ',
	'(',x$dist.method,', ',x$method,')',sep=''),...)
par(old.par)
}

"hclus.table" <-
function(x){

z1<-x$dist.method
z2<-x$method
z3<-seq(length(x$labels)-1,1)
z3<-as.data.frame(cbind(z3,x$merge,x$height))
colnames(z3)<-c('no. clusters','entity','entity','distance')
z<-list(z1,z2,z3)
names(z)<-c('dist.method','method','cluster.table')
return(z)
}

"hist.plots" <-
function(x,var='',by='',save.plot=FALSE,na.rm=TRUE,
	col.hist='blue',col.line='black',las=1,lab=c(5,5,4),...){

old.par<-par(no.readonly=TRUE)

if(!var==''){
	y<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	y<-as.data.frame(y)
	}
else{y<-as.data.frame(x)}

if(by==''){ #histogram w/o groups
	par(mfrow=c(1,1),mar=c(5,5,4,2)) #graphics settings
	for(i in 1:ncol(y)){ #loop thru variables
		par(new=FALSE)
		hist(y[[i]],prob=TRUE,col=col.hist,las=las,lab=lab,
			xaxs='i',yaxs='i',xlab=names(y[i]),
			main=paste('Histogram of',names(y[i]),sep=' '),...)
		par(new=TRUE)
		plot(density(y[[i]],na.rm=na.rm),xaxs='i',yaxs='i',
			col=col.line,las=las,lab=lab,
			axes=FALSE,xlab='',ylab='',main='',...)
		if(save.plot==TRUE){
			dev.print(jpeg,file=paste('hist.',names(y[i]),'.jpg',sep=''),width=800,height=600)
			} #end save
		if(!i==ncol(y)) {readline("Press return for next plot ")}
		} #end loop thru variables
	} #end histogram w/o groups
		
else{ #histograms w/ groups
	n<-by.names(x,by) #create by variable
	y<-cbind(n,y) #bind with selected variables
	groups<-levels(y[,2]) #create vector of group names
	s<-floor(length(groups)^.5) #create multi-figure dimensions
	s<-c(s,ceiling(length(groups)/s)) #create multi-figure dimensions
	for(i in 3:ncol(y)){ #loop thru variables
		par(mfrow=s,mar=c(5,5,4,2)) #graphics settings
		for(j in 1:length(groups)){ #loop thru groups
			z<-y[y[,2]==groups[j],] #select records for group j
			y.min<-min(y[i],na.rm=na.rm)
			y.max<-max(y[i],na.rm=na.rm)
			hist(z[[i]],prob=TRUE,col=col.hist,las=las,lab=lab,
				xaxs='i',yaxs='i',xlab=names(y[i]),main=groups[j],
				xlim=c(y.min,y.max),...)
			par(new=TRUE)
			plot(density(z[[i]],na.rm=na.rm),xaxs='i',yaxs='i',
				col=col.line,las=las,lab=lab,xlim=c(y.min,y.max),
				axes=FALSE,xlab='',ylab='',main='',...)
			if(j==length(groups) & save.plot==TRUE){
				dev.print(jpeg,file=paste('hist.',names(y[i]),'.jpg',sep=''),width=800,height=600)
				} #end save
			} #end loop thru groups
			if(!i==ncol(y)) {readline("Press return for next plot ")}
		} #end loop thru variables
	} #end histogram w/groups
par(old.par)
} #end function

"intrasetcor" <-
function(object){ 

    w <- weights(object)
    lc <- sweep(object$CCA$u, 1, sqrt(w), "*")
    cor(qr.X(object$CCA$QR), lc)
}

"lda.structure" <-
function(x.lda,x,dim=ncol(x.lda),
	digits=3,cutoff=0){

#check for dim limit
if(dim>ncol(x.lda)){
    cat("Only",ncol(x.lda),"axes available\n")
    dim<-length(x.lda)
	}

#calculate structure correlations
z<-cor(x,x.lda[,1:dim])

#print results 
cat("\nStructure Correlations:\n")
z<-round(z,digits=digits)
z[abs(z)<cutoff]<-substring('',1,nchar(z[1,1]))
z<-as.data.frame(z)
return(z)
}

"mantel2" <-
function(xdis,ydis,method='pearson',permutations=1000,strata){

    xdis<-as.dist(xdis)
    ydist<-ydis
    ydis<-as.vector(as.dist(ydis))
    tmp<-cor.test(as.vector(xdis),ydis,method=method)
    statistic<-as.numeric(tmp$estimate)
    variant<-tmp$method
    if(permutations){
        N<-attributes(xdis)$Size
        perm<-rep(0,permutations)
        for(i in 1:permutations){
            take<-permuted.index(N,strata)
            permvec<-as.vector(as.dist(as.matrix(xdis)[take,take]))
            perm[i]<-cor(permvec,ydis,method=method)
			}
        signif<-sum(perm>=statistic)/permutations
	    }
    else{
        signif<-NA
        perm<-NULL
	    }
    res<-list(call=match.call(),method=variant,statistic=statistic, 
        signif=signif,perm=perm,permutations=permutations,xdis=xdis,ydis=ydist)
    if(!missing(strata)){
        res$strata<-deparse(substitute(strata))
        res$stratum.values<-strata
	    }
    class(res)<-'mantel'
    res
}

"mrpp2" <-
function(dat,grouping,permutations=1000,distance='euclidean',
    weight.type=1,strata){

    mrpp.perms<-function(ind,dmat,indls,w){
    	weighted.mean(sapply(indls,function(x) mean(c(dmat[ind==x,ind==x]),na.rm=TRUE)),w=w,na.rm=TRUE)
	    }
    if(inherits(dat,'dist')) 
        dmat<-dat
    else if(is.matrix(dat)&&nrow(dat)==ncol(dat)&&all(dat[lower.tri(dat)]==t(dat)[lower.tri(dat)])){
        dmat<-dat
        attr(dmat,'method')<-'user supplied square matrix'
	    }
    else dmat<-vegdist(dat,method=distance)
    distance<-attr(dmat,'method')
    dmat<-as.matrix(dmat)
    diag(dmat)<-NA
    N<-nrow(dmat)
    ind<-as.numeric(grouping)
    indls<-unique(ind)
    w<-sapply(indls,function(x) sum(ind==x))
    w<-switch(weight.type,w,w-1,w*(w-1)/2)
    del<-mrpp.perms(ind,dmat,indls,w)
    if (missing(strata))
		strata<-NULL
    perms<-sapply(1:permutations,function(x) ind[permuted.index(N,strata=strata)])
    m.ds<-numeric(permutations)
    m.ds<-apply(perms,2,function(x) mrpp.perms(x,dmat,indls,w))
    E.del<-mean(m.ds)
    p<-(1+sum(del>=m.ds))/(permutations+1)
    r2<-1-del/E.del
    
    x<-as.dist(dmat)
	sol<-c(call=match.call())
	grouping<-as.factor(grouping)
	matched<-function(irow,icol,grouping){
	    grouping[irow]==grouping[icol]
		}
	x.vec<-as.vector(x)
	N<-attributes(x)$Size
	irow<-as.vector(as.dist(row(matrix(nrow=N,ncol=N))))
	icol<-as.vector(as.dist(col(matrix(nrow=N,ncol=N))))
	within<-matched(irow,icol,grouping)
	cl.vec<-rep("Between",length(x))
	take<-as.numeric(irow[within])
	cl.vec[within]<-levels(grouping)[grouping[take]]
	cl.vec<-factor(cl.vec,levels=c("Between",levels(grouping)))
    
    out<-list(call=match.call(),delta=del,E.delta=E.del, 
        Pvalue=p,A=r2,distance=distance,weight.type=weight.type, 
        boot.deltas=m.ds,permutations=permutations,class.vec=cl.vec,
        diss.vec=x.vec)
    if (!is.null(strata)){
        out$strata<-deparse(substitute(strata))
        out$stratum.values<-strata
		}
    class(out)<-'mrpp'
    out
}

"multivar.outliers" <-
function(x,method,var='',cor.method='pearson',
	outfile='',sd.limit=3,alpha=.001,plot=TRUE,save.plot=FALSE,
	use='complete.obs',na.rm=TRUE,col.hist='blue',col.line='black',
	col.point='blue',las=1,lab=c(5,5,4),...){

library(vegan)
library(MASS)

old.par<-par(no.readonly=TRUE)

if(!var==''){
	y<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	}
else{y<-x}

#paired histogram function for avedist and sd(avedist)
paired.hist<-function(file=d5){
	par(mfrow=c(2,1),mar=c(5,5,4,2)) #graphics settings
	hist(file[[1]],prob=TRUE,col=col.hist,las=las,lab=lab,
		xaxs='i',yaxs='i',xlab='average distance',
		main=paste('Histogram of',method,'distance',sep=' '),...)
	par(new=TRUE)
	plot(density(file[[1]],na.rm=na.rm),
		xaxs='i',yaxs='i',col=col.line,las=las,lab=lab,
		axes=FALSE,xlab='',ylab='',main='',...)
	hist(file[[2]],prob=TRUE,col=col.hist,las=las,lab=lab,
		xaxs='i',yaxs='i',xlab='standard deviation of average distance',
		main='',...)
	par(new=TRUE)
	plot(density(file[[2]],na.rm=na.rm),
		xaxs='i',yaxs='i',col=col.line,las=las,lab=lab,
		axes=FALSE,xlab='',ylab='',main='',...)
	if(save.plot==TRUE){ #save plot to file
		dev.print(jpeg,file=paste('dhist.',method,'.jpg',sep=''),width=800,height=600)
		} #end save	plot
	} #end paired histogram function

#ecological distance calculations
if(method=='mahalanobis'){ #compute mahalanobis distance
	stopifnot(mahalanobis(y,0,diag(ncol(y)))==rowSums(y*y)) #Here,D^2=usual squared Euclidean distances
#	y.pca<-prcomp(y,center=TRUE,scale=TRUE)
#	Sy<-cov.rob(y.pca$x)
	Sy<-cov.rob(y)
	z<-mahalanobis(y,Sy$center,Sy$cov)
	n<-nrow(x)
	p<-ncol(x)
	critical.chi<-qchisq(1-alpha,p) # asymptotic chi-square method	 
	if(plot==TRUE){
		par(mfrow=c(2,1),mar=c(5,5,4,2)) #graphics settings
		plot(density(z),col=col.line,las=las,lab=lab,
			main=paste('Squared Mahalanobis distances', 
			'(n=',nrow(y),'p=',ncol(y),')',sep=' '),...)
		rug(z)
		abline(v=critical.chi,col='red',lty=2)
		qqplot(qchisq(ppoints(nrow(y)),df=ncol(y)),z,
			col=col.point,las=las,lab=lab,
			xlab='Theoretical Quantiles',ylab='Sample Quantiles',
			main=expression("Q-Q plot of Mahalanobis" * ~D^2 *
	        " vs. quantiles of" * ~ chi[3]^2),...)
	    abline(0,1,col='red',lty=2)
	    }
    z<-as.data.frame(z)
    names(z)<-c('mahalanobis distance')
	t1<-z>=critical.chi | z<=-sd.limit #identify extreme values
	t2<-as.matrix(t1)
	t3<-z[apply(t2,1,FUN='any',na.rm=na.rm),,drop=FALSE] #select rows with extremes
	z<-round(t3,3) #round results to 3 decimal places
	if(!outfile==''){ #save outfile
		write.table(z,file=paste(outfile,'.csv',sep=''),quote=FALSE,sep=',')
		} #end save outfile
 	} #end mahalanobis distance

else{ #all other distances 
	if(method=='correlation'){ #compute correlation distance
		d1<-as.matrix((1-cor(t(y),method=cor.method,use=use))/2) #compute distance
		} #end correlation distance
	else{ #compute other distance methods
		d1<-as.matrix(vegdist(y,method=method,na.rm=na.rm,...)) #vegdist function
		} #end other distance methods
	diag(d1)<-NA
	d2<-apply(d1,2,mean,na.rm=na.rm)
	d3<-(d2-mean(d2))/sd(d2)
	d4<-c('avedist','sd')
	d5<-cbind(d2,d3)
	colnames(d5)<-d4
	d5<-as.data.frame(d5)
	if(plot==TRUE){ #plot paired histograms
		paired.hist(d5)
		} #end plot paired histograms
	t1<-d3>=sd.limit | d3<=-sd.limit #identify extreme values
	t2<-as.matrix(t1)
	t3<-d5[apply(t2,1,FUN='any',na.rm=na.rm),,drop=FALSE] #select rows with extremes
	z<-round(t3,3) #round results to 3 decimal places
	if(!outfile==''){ #save outfile
		write.table(z,file=paste(outfile,'.csv',sep=''),quote=FALSE,sep=',')
		} #end save outfile
	} #end all other distances
	
par(old.par)
return(z)
} #end function	

"nhclus.scree" <-
function(x,max.k,...){

library(cluster)
old.par<-par(no.readonly=TRUE)

#nonhierarchical clustering
s<-rep(0,max.k)
d<-rep(0,max.k)
for(i in 2:max.k){
	temp<-pam(x,k=i,...)
	s[i]<-temp$silinfo$avg.width #get result (ave silhouette width)
	d[i]<-temp$objective[2] #get result (min sum of diss to medoid)
	}

#merge results into data.frame
s<-s[-1]
d<-d[-1]
n<-c(2:max.k)
y<-round(as.data.frame(cbind(n,d,s)),3)
names(y)<-c('no. clusters','sum within-clus diss','ave silhouette width')

#create scree plot
par(mar=c(5.1,4.5,4.1,4.5))
plot(y[,1],y[,2],type='o',lwd=1.5,pch=19,col='blue',
	ylab='Within-cluster Dissimilarity',xlab='Number of Clusters')
par(new=TRUE)
plot(y[,1],y[,3],type='o',lty=2,lwd=1.5,pch=19,col='red',
	xlab='',ylab='',yaxt='n')
	sil.tick<-pretty(range(y[,3]))
	axis(side=4,at=sil.tick,srt=90)
	mtext('Average Silhouette Width',side=4,line=3)
legend(x='topright',inset=c(.05,.05),legend=c('Within-cluster dissimilarity','Average silhoutte width'),col=c('blue','red'),lty=c(1,3),lwd=c(1.5,1.5))
title(main='Scree Plot of K-means Clustering')

par(old.par)
return(y)
}

"nmds.monte" <-
function(x,k,distance='bray',trymax=50,autotransform=FALSE,
	trace=0,zerodist='add',perm=100,col.hist='blue',col.line='red',
	lty=2,las=1,lab=c(5,5,4),...){

library(vegan)
library(MASS)

z<-metaMDS(comm=x,k=k,distance=distance,trymax=trymax,autotransform=autotransform,trace=trace,...) #nmds analysis
z.stress<-z$stress #get stress
y.stress<-rep(0,perm)
for(i in 1:perm){
	y<-apply(comm,2,sample) #permute data matrix
	y<-metaMDS(comm=y,k=k,distance=distance,trymax=trymax,autotransform=autotransform,trace=trace,...) #nmds analysis
	y.stress[i]<-y$stress #get stress
	}
n<-sum(y.stress<=z.stress) #compute number of random runs with stress < observed
p.value<-(1+n)/(1+perm) #compute p-value

xmin<-min(z.stress,min(y.stress))
xmax<-max(z.stress,max(y.stress))
hist(y.stress,col=col.hist,las=las,lab=lab,
	xaxs='i',yaxs='i',xlim=c(xmin,xmax),xlab='Stress',
	main=paste('Random Permutation Distribution of Stress for',k,'Dimensions',sep=' '),...)
abline(v=z.stress,col=col.line,lty=lty,lwd=2,...)

cat('Randomization Test of Stress:\n')
cat('Permutation stress values:\n')
print(y.stress)
z<-rbind('Observed stress'=z.stress,'P-value'=p.value)
return(z)
}

"nmds.scree" <-
function(x,distance='bray',k=6,trymax=50,
	autotransform=FALSE,trace=0,...){

library(vegan)
library(MASS)

old.par<-par(no.readonly=TRUE)

nmds.stress<-rep(0,k)
nmds.dim<-c(1:k)
	for(i in 1:k){
	y<-metaMDS(x,distance=distance,k=i,trymax=trymax,
		autotransform=autotransform,trace=trace,...)
	nmds.stress[i]<-y$stress
	}
plot(nmds.dim,nmds.stress,type='o',pch=19,col='blue',
	ylab='Stress',xlab='Ordination Axis',
	main='Scree Plot of Stress vs. Dimension',...)

par(old.par)
}

"norm.test" <-
function(x,groups='',var='',method='ad',...){

library(nortest)

if(!var==''){
	y<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	y<-as.data.frame(y)
	}
else{y<-as.data.frame(x)}

if(!groups==''){
	y.resid<-matrix(0,nrow(y),ncol(y))
	grp<-as.factor(groups)
	for(i in 1:ncol(y)){ #loop thru variables
		y.aov<-aov(y[,i]~grp)
		y.resid[,i]<-y.aov$residuals
		}
	y.resid<-as.data.frame(y.resid)
	rownames(y.resid)<-rownames(y)
	colnames(y.resid)<-colnames(y)
	y<-y.resid
	}
	
#create summary table
z<-matrix(0,ncol(y),2) #create blank matrix
for(i in 1:ncol(y)){ #loop thru variables
	if(method=='ad'){
		temp<-ad.test(y[,i],...)
		}
	else if(method=='sf'){
		temp<-sf.test(y[,i],...)
		}		
	if(method=='cvm'){
		temp<-cvm.test(y[,i],...)
		}
	else if(method=='lillie'){
		temp<-lillie.test(y[,i],...)
		}
	else if(method=='pearson'){
		temp<-pearson.test(y[,i],...)
		}
	z[i,1]<-round(temp$statistic,3)
	z[i,2]<-temp$p.value
	}
rownames(z)<-colnames(y)
if(method=='ad'){
	colnames(z)<-c('Anderson-Darling A','p-value')
	cat('Anderson-Darling Test of Normality:\n')
	}
if(method=='sf'){
	colnames(z)<-c('Shapiro-Francia','p-value')
	cat('Shapiro-Francia Test of Normality:\n')
	}
if(method=='cvm'){
	colnames(z)<-c('Cramer-von Mises W','p-value')
	cat('Cramer-van Mises Test of Normality:\n')
	}
if(method=='lillie'){
	colnames(z)<-c('Lilliefors D','p-value')
	cat('Lilliefors (Kolmogorov-Smirnov) Test of Normality:\n')
	}
if(method=='pearson'){
	colnames(z)<-c('Pearson chi-square','p-value')
	cat('Pearson chi-square Test of Normality:\n')
	}
z<-as.data.frame(z)
z[,2]<-format.pval(z[,2],digits=3,eps=.001)		
return(z)
}

"ordi.monte" <-
function(x,ord,dim=length(x),perm=1000,center=TRUE,
	scale=TRUE,digits=3,plot=TRUE,col.hist='blue',col.line='red',
	lty=2,las=1,lab=c(5,5,4),...){

p<-length(x)
if(dim>p){
    cat("Only",p,"axes available\n")
    dim<-p
	}

if(ord=='pca'){
	z<-prcomp(x,center=center,scale=scale) #prcomp analysis
	z.eig<-z$sdev[1:dim]^2 #compute eigenvalues
	z.teig<-t(z.eig) #transpose eigenvalues
	z.teig<-t(matrix(z.teig,length(z.teig),perm)) #fill matrix with eigenvalues
	write('',file='y.csv') #empty outfile if it exists
	for(i in 1:perm){
		y<-apply(x,2,sample) #permute data matrix
		y<-prcomp(y,center=center,scale=scale) #prcomp on permuted matrix
		y<-y$sdev[1:dim]^2 #compute eigenvalues
		y<-as.data.frame(t(y)) #coerce to data frame and transpose
		write.table(y,file='y.csv',sep=',',append=TRUE,row.names=FALSE,col.names=FALSE)
		}
	y<-read.table('y.csv',sep=',',header=FALSE) #read in permutation results
	p.value<-apply(y>z.teig,2,sum) #compute proportion of random distribution smaller than observed
	p.value<-p.value/perm #compute p-value
	names<-colnames(z$rotation[,1:dim]) #add 'PC#' names
	}

else if(ord=='ca'){
	library(vegan)
	z<-cca(x) #correspondence analysis
	z.eig<-z$CA$eig[1:dim] #get eigenvalues
	z.teig<-t(z.eig) #transpose eigenvalues
	z.teig<-t(matrix(z.teig,length(z.teig),perm)) #fill matrix with eigenvalues
	write('',file='y.csv') #empty outfile if it exists
	for(i in 1:perm){
		y<-apply(x,2,sample) #permute data matrix
		y<-cca(y) #CA on permuted matrix
		y<-y$CA$eig[1:dim] #get eigenvalues
		y<-as.data.frame(t(y)) #coerce to data frame and transpose
		write.table(y,file='y.csv',sep=',',append=TRUE,row.names=FALSE,col.names=FALSE)
		}
	y<-read.table('y.csv',sep=',',header=FALSE) #read in permutation results
	p.value<-apply(y>z.teig,2,sum) #compute proportion of random distribution smaller than observed
	p.value<-p.value/perm #compute p-value
	names<-names(z$CA$eig[1:dim]) #add 'CA#' names
	}

else if(ord=='dca'){
	library(vegan)
	if(dim>4){
    cat("Only",4,"axes available\n")
    dim<-4
	}
	z<-decorana(x,...) #detrended correspondence analysis
	z.eig<-z$evals[1:dim] #get eigenvalues
	z.teig<-t(z.eig) #transpose eigenvalues
	z.teig<-t(matrix(z.teig,length(z.teig),perm)) #fill matrix with eigenvalues
	write('',file='y.csv') #empty outfile if it exists
	for(i in 1:perm){
		y<-apply(x,2,sample) #permute data matrix
		y<-decorana(y,...) #DCA on permuted matrix
		y<-y$evals[1:dim] #get eigenvalues
		y<-as.data.frame(t(y)) #coerce to data frame and transpose
		write.table(y,file='y.csv',sep=',',append=TRUE,row.names=FALSE,col.names=FALSE)
		}
	y<-read.table('y.csv',sep=',',header=FALSE) #read in permutation results
	p.value<-apply(y>z.teig,2,sum) #compute proportion of random distribution smaller than observed
	p.value<-p.value/perm #compute p-value
	names<-names(z$eval[1:dim]) #add 'CA#' names
	}

if(plot==TRUE){
	for(i in 1:dim){
		xmin<-min(min(y[[i]],z.eig[i]))
		xmax<-max(max(y[[i]],z.eig[i]))
		hist(y[[i]],col=col.hist,las=las,lab=lab,
			xaxs='i',yaxs='i',xlim=c(xmin,xmax),
			xlab='Eigenvalue',
			main=paste('Random Permutation Distribution of Eigenvalues for',names[i],sep=' '),...)
		abline(v=z.eig[i],col=col.line,lty=lty,lwd=2,...)
	    readline("Press return for next plot ")
		}
	}	

cat('Randomization Test of Eigenvalues:\n')
z<-rbind('Eigenvalue'=z.eig,'P-value'=p.value)
colnames(z)<-names
z<-round(z,digits=digits)
return(z)
}

"ordi.overlay" <-
function(x.ord,x,var='',fit=TRUE,choices=c(1,2),expand=5,
	alpha=.95,pch=19,...){

library(fields)

if(!var==''){
	y<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	}
else{y<-x}

old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))
scores<-scores(x.ord,display='site',choices=choices)
names<-colnames(scores)

for(i in 1:length(y)){
	if(fit==TRUE){
		set.panel(2,2,relax=TRUE)
	
		#top left plot
		out<-qsreg(scores[,2],y[,i],alpha=alpha)
	    plot(out$x,out$y,xlab=names[2],ylab=names(y[i]),pch=pch)
	    orderx<-order(out$x)
	    temp<-out$fitted.values[,c(out$ind.cv,out$ind.cv.ps)]
	    matlines(out$x[orderx],temp[orderx,],lty=1,col=c('blue','red'))	
		
		#top right plot
		plot(scores[,1],scores[,2],col='blue',pch=pch,
			xlab=names[1],ylab=names[2],
			cex=expand*y[,i]/max(y[,i]),main=names(y[i]),...)
	
		#bottom left plot
		frame()
		
		#bottom right plot
		out<-qsreg(scores[,1],y[,i],alpha=alpha)
	    plot(out$x,out$y,xlab=names[1],ylab=names(y[i]),pch=pch)
	    orderx<-order(out$x)
	    temp<-out$fitted.values[,c(out$ind.cv,out$ind.cv.ps)]
	    matlines(out$x[orderx],temp[orderx,],lty=1,col=c('blue','red'))	
	
		readline("Press return for next plot ")    
		}
	else{
		plot(scores[,1],scores[,2],col='blue',pch=pch,
			xlab=names[1],ylab=names[2],
			cex=expand*y[,i]/max(y[,i]),main=names(y[i]))
		readline("Press return for next plot ")
		}
	}
}

"ordi.scree" <-
function(x,ord,...){

library(vegan)
library(MASS)

old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))

if(ord=='pca'){
	if(class(x)=='prcomp'||class(x)=='princomp'){
		eig<-x$sdev^2
		if(!x$scale==''){
			inertia<-'correlations'
			}
		else{
			inertia<-'variance'
			}
		}
	else if(any(class(x)=='rda')){ #from rda()
		eig<-x$CA$eig
		inertia==x$inertia
		}
	#calculate broken stick values for PCA
	if(inertia=='correlations'){
		p<-length(eig)
		y<-rep(0,p)
		for(i in 1:p) y[i]<-1/i
		for(i in 1:p) y[i]<-sum(y[i:p]) #broken stick values
		}
	}
else if(ord=='ca'){ #from cca()
	eig<-x$CA$eig
	}
else if(ord=='mds'){ #from cmdscale()
	eig<-x$eig
	}			
else if(ord=='rda'||ord=='cca'||ord=='cmds'){ #from rda(),cca(),or capscale()
	eig<-x$CCA$eig
	}
	
#create scree plot
if(ord=='pca'&&inertia=="correlations"){ #PCA with broken stick
	plot(eig,type='o',lwd=1.5,pch=19,col='blue',ylim=c(0,max(max(eig),max(y))),
		ylab='Eigenvalue',xlab='Ordination Axis',
		main='Scree Plot of Unconstrained Eigenvalues',...)
	par(new=TRUE)
	plot(y,type='o',lty=3,lwd=1.5,pch=19,col='green',ylim=c(0,max(max(eig),max(y))),
		xlab='',ylab='',main='')
	abline(h=sum(eig)/p,col='red')
	legend(x='topright',inset=c(.1,.1),legend=c('Observed','Broken-stick','Latent root'),col=c('blue','green','red'),lty=c(1,3,1))
	}
else{ #all others
	plot(eig,type='o',lwd=1.5,pch=19,col='blue',ylim=c(0,max(eig)),
		ylab='Eigenvalue',xlab='Ordination Axis',...)
		if(ord=='pca'||ord=='ca'||ord=='mds'){
			title(main='Scree Plot of Unconstrained Eigenvalues')
			}
		else if(ord=='rda'||ord=='cca'||ord=='cmds'){
			title(main='Scree Plot of Constrained Eigenvalues')
			}
	}
#create cumulative scree plot	
readline("Press return for next plot ")	
plot(cumsum(eig/sum(eig)),type='o',pch=19,col='blue',
	ylab='Cumulative Proportion',xlab='Ordination Axis',...)
	if(ord=='pca'||ord=='ca'||ord=='mds'){
		title(main='Cumulative Scree Plot of Unconstrained Eigenvalues')
		}
	else if(ord=='rda'||ord=='cca'||ord=='cmds'){
		title(main='Cumulative Scree Plot of Constrained Eigenvalues')
		}	
}

"pca.communality" <-
function(x.pca,x,dim=length(x.pca$sdev),
	digits=3){
#don't know why communality for all dimensions doesn't always
#equal one, perhaps singularity issue.

#check for dim limit
if(dim>length(x.pca$sdev)){
    cat("Only",length(x.pca$sdev),"axes available\n")
    dim<-length(x.pca$sdev)
	}

#calculate final communality estimates
z<-cor(x,x.pca$x[,1:dim])
z<-z^2
c<-apply(z,1,sum)

#print results 
cat("\nFinal Communalities -",dim,'dimensions:\n',sep=' ')
c<-round(c,digits=digits)
c<-as.data.frame(c)
return(c)
}

"pca.eigenval" <-
function(x.pca,dim=length(x.pca$sdev),digits=7){

#check for dim limit
if(dim>length(x.pca$sdev)){
    cat("Only",length(x.pca$sdev),"axes available\n")
    dim<-length(x.pca$sdev)
	}

#calculate some variables
names<-colnames(x.pca$rotation[,1:dim])
var<-x.pca$sdev^2
trace<-sum(var)
prop.var<-var/trace

#broken-stick distribution
p<-length(x.pca$sdev)
y<-rep(0,p)
for(i in 1:p) y[i]<-1/i
for(i in 1:p) y[i]<-sum(y[i:p])
y<-y[1:dim]

#print results    
cat('Importance of components:\n')
z<-rbind('Variance(eigenvalue)'=var[1:dim],
	'Proportion of Variance'=prop.var[1:dim],
	'Cumulative Proportion'=cumsum(prop.var[1:dim]),
	'Broken-stick value'=y)
colnames(z)<-names
z<-round(z,digits=digits)
return(z)
}

"pca.eigenvec" <-
function(x.pca,dim=length(x.pca$sdev),
	digits=7,cutoff=0){

#check for dim limit	
if(dim>ncol(x.pca$rotation)){
    cat("Only",ncol(x.pca$rotation),"axes available\n")
    dim<-ncol(x.pca$rotation)
    }

#print results    
cat("\nEigenvectors:\n")
z<-format(round(x.pca$rotation[,1:dim],digits=digits))
z[abs(x.pca$rotation[,1:dim])<cutoff]<-substring('',1,nchar(z[1,1]))
z<-as.data.frame(z)
return(z)
}

"pca.structure" <-
function(x.pca,x,dim=length(x.pca$sdev),
	digits=3,cutoff=0){

#check for dim limit
if(dim>length(x.pca$sdev)){
    cat("Only",length(x.pca$sdev),"axes available\n")
    dim<-length(x.pca$sdev)
	}

#calculate structure correlations
z<-cor(x,x.pca$x[,1:dim])

#print results 
cat("\nStructure Correlations:\n")
z<-round(z,digits=digits)
z[abs(z)<cutoff]<-substring('',1,nchar(z[1,1]))
z<-as.data.frame(z)
return(z)
}

"plot.anosim" <-
function(x,title1='ANOSIM (within- vs between-group rank dissimilarities)',
	title2='ANOSIM (observed vs expected R)',col='blue',...){

old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))

#boxplot between and within-group dissimilarities
boxplot(x$dis.rank~x$class.vec,notch=TRUE,varwidth=TRUE,col=col,
	ylab='Rank Dissimilarity',...)
title(title1)
if(x$permutations){
	pval<-format.pval(x$signif,eps=1/x$permutations)
	}
else{
	pval<-'not assessed'
	}
mtext(paste("R = ",round(x$statistic,3),', ','P = ',pval),3)

readline("Press return for next plot ")

#histogram of permutation distribution
r<-x$statistic
E.r<-mean(x$perm)
perm<-x$perm
hist(perm,col=col,xlim=range(c(r,perm)),xlab='R',main='')
abline(v=r,col='red',lwd=2,lty=3)
title(title2)
if(x$permutations){
	pval<-format.pval(x$signif,eps=1/x$permutations)
	}
else{
	pval<-'not assessed'
	}
mtext(paste('Observed R = ',round(x$statistic,3),',',
	'Expected R = ',round(E.r,3),',','P = ',pval),3)

}

"plot.mantel" <-
function(x,title1='MANTEL Scatterplot)',
	title2='MANTEL (observed vs expected R)',col='blue',...){

old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))

#scatterplot of dissimilarity matrices
plot(x$xdis,x$ydis,col=col,xlab='X-matrix dissimilarities',
	ylab='Y-matrix dissimilarities',...)
title(title1)
if(x$permutations){
	pval<-format.pval(x$signif,eps=1/x$permutations)
	}
else{
	pval<-'not assessed'
	}
mtext(paste("R = ",round(x$statistic,3),', ','P = ',pval),3)

readline("Press return for next plot ")

#histogram of permutation distribution
r<-x$statistic
E.r<-mean(x$perm)
perm<-x$perm
hist(perm,col=col,xlim=range(c(r,perm)),xlab='r',main='')
abline(v=r,col='red',lwd=2,lty=3)
title(title2)
if(x$permutations){
	pval<-format.pval(x$signif,eps=1/x$permutations)
	}
else{
	pval<-'not assessed'
	}
mtext(paste('Observed r = ',round(x$statistic,3),',',
	'Expected r = ',round(E.r,3),',','P = ',pval),3)

}

"plot.mrpp" <-
function(x,title1='MRPP (within- vs between-group dissimilarities)',
	title2='MRPP (observed vs expected delta)',col='blue',...){

old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))

#boxplot of between and within-group dissimilarities
boxplot(x$diss.vec~x$class.vec,notch=TRUE,varwidth=TRUE,col=col,
	ylab='Dissimilarity',...)
title(title1)
if(x$permutations){
	pval<-format.pval(x$Pvalue,eps=1/x$permutations)
	}
else{
	pval<-'not assessed'
	}
mtext(paste('Delta = ',round(x$delta,3),',','P = ',pval),3)

readline("Press return for next plot ")

#histogram of permutation distribution
d<-x$delta
bd<-x$boot.deltas
hist(bd,col=col,xlim=range(c(d,bd)),xlab='Delta',main='')
abline(v=d,col='red',lwd=2,lty=3)
title(title2)
if(x$permutations){
	pval<-format.pval(x$Pvalue,eps=1/x$permutations)
	}
else{
	pval<-'not assessed'
	}
mtext(paste('Observed delta = ',round(x$delta,3),',',
	'Expected delta = ',round(x$E.delta,3),',','P = ',pval),3)

}

"qqnorm.plots" <-
function(x,var='',by='',save.plot=FALSE,
	na.rm=TRUE,col.point='blue',col.line='red',las=1,...){

old.par<-par(no.readonly=TRUE)
on.exit(par(old.par))

if(!var==''){
	y<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	y<-as.data.frame(y)
	}
else{y<-as.data.frame(x)}

if(by==''){ #QQnorm plots w/o groups
	for(i in 1:ncol(y)){ #loop thru variables
		qqnorm(y[,i],datax=TRUE,col=col.point,las=las,
			main=paste('Normal Q-Q Plot of',names(y[i]),sep=' '),...)
		y.IQR<-IQR(y[,i],na.rm=na.rm)
		if(y.IQR>0){ #add qqline
			qqline(y[,i],datax=TRUE,col=col.line,...)
			} #end qqline
		if(save.plot==TRUE){ #save plot
			dev.print(jpeg,file=paste('qqnorm.',names(y[i]),'.jpg',sep=''),width=800,height=600)
			} #end save plot
		if(!i==ncol(y)) {readline("Press return for next plot ")}
		} #end loop thru variables
	} #end qqnorm w/o groups
	
else{ #QQnorm plots w/ groups
	n<-by.names(x,by) #create by variable
	y<-cbind(n,y) #bind with selected variables
	groups<-levels(y[,2]) #create vector of group names
	s<-floor(length(groups)^.5) #create multi-figure dimension
	s<-c(s,ceiling(length(groups)/s)) #create multi-figure dimensions
	for(i in 3:ncol(y)){ #loop thru variables
		par(mfrow=s,mar=c(5,5,4,2))	#graphics settings	
		for(j in 1:length(groups)){ #loop thru groups
			z<-y[y[,2]==groups[j],] #select records for group j
			qqnorm(z[,i],datax=TRUE,col=col.point,las=las,
				main=groups[j],
				ylab=paste('Sample quantiles of',names(z[i]),sep=' '),...)
			z.IQR<-IQR(z[,i],na.rm=na.rm)
			if(z.IQR>0){ #add qqline
				qqline(z[,i],datax=TRUE,col=col.line,...)
				} #end qqline
			if(j==length(groups) & save.plot==TRUE){
				dev.print(jpeg,file=paste('qqnorm.',names(z[i]),'.jpg',sep=''),width=800,height=600)
				} #end save
			} #end loop thru groups
			if(!i==ncol(y)) {readline("Press return for next plot ")}
		} #end loop thru variables
	} #end qqnorm w/ groups
} #end function

"ran.split" <-
function(x,grouping='',prop=.5){

if(!grouping==''){
	y<-cbind(grouping,x)
	}
else{
	y<-x #groups assumed to be in first column
	}
	
N<-nrow(y)
g<-runif(N)
g[g<prop]<-0
g[g>=prop]<-1
y<-cbind(g,y)
y<-split(y,g)
y1<-y[[1]]
y2<-y[[2]]
calibrate<<-y1[,-c(1,2)]
validate<<-y2[,-c(1,2)]
grp.cal<<-y1[,2]
grp.val<<-y2[,2]

z1<-c(nrow(y1),nrow(y2))
z2<-round(c(nrow(y1)/N,nrow(y2)/N),2)
z1<-cbind(z1,z2)
colnames(z1)<-c('Number of samples','Proportion')
rownames(z1)<-c('Calibration set','Validation set')
z2<-table(grp.cal)
z3<-table(grp.val)
z<-list(z1,z2,z3)
names(z)<-c('Random Subset Summary:','Calibration Table','Validation Table')
return(z)
}

"redun.plot" <-
function(x,var='',...){

old.par<-par(no.readonly=TRUE)

if(!var==''){
	z.obs<-as.matrix(cor(x))
	diag(z.obs)<-NA
	z.obs<-as.data.frame(z.obs)

	for(i in 1:ncol(x)){
		y<-apply(x,2,sample) #permute data matrix
		}
	z.ran<-as.matrix(cor(y))
	diag(z.ran)<-NA
	z.ran<as.data.frame(z.ran)

	y.obs<-subset(z.obs,select=eval(parse(text=var))) #select variables to plot
	y.ran<-subset(z.ran,select=eval(parse(text=var))) #select variables to plot

	for(i in 1:ncol(y.obs)){ #loop thru variables
		par(new=FALSE)
		plot(sort(y.obs[,i]),type='l',lwd=2,col='blue',
			xaxs='i',yaxs='i',ylim=c(-1,1),
			xlab='Rank order of pairwise correlations',
			ylab='Correlation',
			main=paste('Redundancy of ',names(y.obs[i]),' vs. random data',sep=''),...)
		par(new=TRUE)
		plot(sort(y.ran[,i]),type='l',lwd=2,col='green',
			xaxs='i',yaxs='i',ylim=c(-1,1),
			xlab='',
			ylab='',
			main='',...)
		abline(h=0,col='red',lty=3,lwd=1,...)
		legend(x='bottomright',inset=c(.1,.1),legend=c('Actual','Random'),
			col=c('blue','green'),lty=c(1,1),lwd=c(2,2))
		if(!i==ncol(y)) {readline("Press return for next plot ")}
		} #end loop thru variables
	}

else{
	z.obs<-sort(as.vector(as.dist(cor(x))))
	for(i in 1:ncol(x)){
		y<-apply(x,2,sample) #permute data matrix
		}
	z.ran<-sort(as.vector(as.dist(cor(y))))
	plot(z.obs,type='l',lwd=2,col='blue',
		xaxs='i',yaxs='i',ylim=c(-1,1),
		xlab='Rank order of pairwise correlations',
		ylab='Correlation',
		main='Redundancy of actual vs. random data',...)
	par(new=TRUE)
	plot(z.ran,type='l',lwd=2,col='green',
		xaxs='i',yaxs='i',ylim=c(-1,1),
		xlab='',
		ylab='',
		main='',...)
	abline(h=0,col='red',lty=3,lwd=1,...)
	legend(x='bottomright',inset=c(.1,.1),legend=c('Actual','Random'),
		col=c('blue','green'),lty=c(1,1),lwd=c(2,2))
	}

par(old.par)
}

"replace.missing" <-
function(x,var='',method='median',outfile=''){

if(!var==''){
	y1<-subset(x,select=-eval(parse(text=var))) #select all other variables
	y2<-subset(x,select=eval(parse(text=var))) #select all specified variables
	}
else{
	y2<-x #original variables
	}

if(method=='mean'){ #for method=mean
	for(i in names(y2)){ #loop through selected variables
		t<-round(mean(y2[[i]],na.rm=TRUE),3) #compute mean for each variable
		y2[is.na(y2[[i]]),i]<-t #assign mean value to missing value
		} 
	}
if(method=='median'){ #for method=median
	for(i in names(y2)){ #loop through selected variables
		t<-median(y2[[i]],na.rm=TRUE) #compute median for each variable
		y2[is.na(y2[[i]]),i]<-t #assign median value to missing value
		} 
	}

if(!var==''){
	z<-cbind(y1,y2) #combine deselected and (modified) selected variables
	}
else{z<-y2}
if(!outfile==''){ #save outfile
	write.table(z,file=paste(outfile,'.csv',sep=''),quote=FALSE,sep=',')
	} #end save outfile
return(z)
} #end function

"sum.stats" <-
function(x,var='',by='',margin='column',...){

if(!var==''){
	y<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	}
else{y<-x}

variable<-colnames(y)
sample<-rownames(y)

#statistical functions
nobs<<-function(x,na.rm) length(x)
cv<<-function(x,na.rm) sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)*100 
xeros<<-function(x,na.rm) sum(x==0,na.rm=TRUE)
pct.xeros<<-function(x,na.rm) sum(x==0,na.rm=TRUE)/length(x)*100
nobs.missing<<-function(x,na.rm) sum(is.na(x))
pct.missing<<-function(x,na.rm) sum(is.na(x))/length(x)*100 
se<<-function(x,na.rm) sd(x,na.rm=TRUE)/sqrt(length(x)-sum(is.na(x)))
se.ratio<<-function(x,na.rm) se(x)/mean(x,na.rm=TRUE)*100
richness<<-function(x,na.rm) nobs(x)-xeros(x)-nobs.missing(x)
sh.diversity<<-function(x,na.rm) -sum(((x)/sum(x,na.rm=TRUE))*log(((x)/sum(x,na.rm=TRUE))),na.rm=TRUE)
sh.evenness<<-function(x,na.rm) sh.diversity(x)/log(richness(x))
si.diversity<<-function(x,na.rm) 1-sum(((x)/sum(x,na.rm=TRUE))*((x)/sum(x,na.rm=TRUE)),na.rm=TRUE)
si.evenness<<-function(x,na.rm) si.diversity(x)/(1-(1/richness(x)))

if(by==''){ #summary table w/o groups
	if(margin=='column'){ #summary table by columns
		z1<-data.frame(apply(y,2,function(x){ #calculate stats
			z1<-c(nobs(x),min(x,na.rm=TRUE),max(x,na.rm=TRUE),
		    mean(x,na.rm=TRUE),median(x,na.rm=TRUE),sum(x,na.rm=TRUE),
		    sd(x,na.rm=TRUE),cv(x),xeros(x),pct.xeros(x),nobs.missing(x),
			pct.missing(x),se(x),se.ratio(x),richness(x),sh.diversity(x),
			sh.evenness(x),si.diversity(x),si.evenness(x))
		    names(z1)<-c('nobs','min','max','mean',
		    'median','sum','sd','cv','xeros','pct.xeros',
		    'nobs.missing','pct.missing','se','se.ratio',
		    'richness','sh.diversity','sh.evenness',
		    'si.diversity','si.evenness') #create col names
			z1<-round(z1,3) #round elements to 3 decimal places
			}))
		z2<-data.frame(t(apply(z1,1,function(x){ #calculate stats on stats
			z2<-c(nobs(x),min(x,na.rm=TRUE),max(x,na.rm=TRUE),
			mean(x,na.rm=TRUE),median(x,na.rm=TRUE),sd(x,na.rm=TRUE),cv(x))
		    names(z2)<-c('nobs','min','max','mean',
		    'median','sd','cv') #create row names
			z2<-round(z2,3) #round elements to 3 decimal places
			})))
		z<-list(z1,z2) #create list with col stats and sum stats
		names(z)<-c('Column.Summary','Table.Summary')
		} #end summary table by columns

	else{ #summary table by rows
		z1<-data.frame(t(apply(y,1,function(x){ #calculate stats
			z1<-c(nobs(x),min(x,na.rm=TRUE),max(x,na.rm=TRUE),
		    mean(x,na.rm=TRUE),median(x,na.rm=TRUE),sum(x,na.rm=TRUE),
		    sd(x,na.rm=TRUE),cv(x),xeros(x),pct.xeros(x),nobs.missing(x),
			pct.missing(x),se(x),se.ratio(x),richness(x),sh.diversity(x),
			sh.evenness(x),si.diversity(x),si.evenness(x))
		    names(z1)<-c('nobs','min','max','mean',
		    'median','sum','sd','cv','xeros','pct.xeros',
		    'nobs.missing','pct.missing','se','se.ratio',
		    'richness','sh.diversity','sh.evenness',
		    'si.diversity','si.evenness') #create col names
			z1<-round(z1,3) #round elements to 3 decimal places
			})))
		z2<-data.frame(apply(z1,2,function(x){ #calculate stats on stats
			z2<-c(nobs(x),min(x,na.rm=TRUE),max(x,na.rm=TRUE),
			mean(x,na.rm=TRUE),median(x,na.rm=TRUE),sd(x,na.rm=TRUE),cv(x))
		    names(z2)<-c('nobs','min','max','mean',
		    'median','sd','cv') #create row names
			z2<-round(z2,3) #round elements to 3 decimal places
			}))
		z<-list(z1,z2) #create list with row stats and sum stats
		names(z)<-c('Row.Summary','Table.Summary')
		} #end summary table by rows
	} #end summary table w/o groups

else{ #summary table w/ groups
#	write('',file=paste(outfile,'.csv',sep='')) #empty outfile if it exists
	fns<-c('nobs','min','max','mean',
	    'median','sum','sd','cv','xeros','pct.xeros',
	    'nobs.missing','pct.missing','se','se.ratio',
	    'richness','sh.diversity','sh.evenness',
	    'si.diversity','si.evenness') #create names vector
	n<-by.names(x,by) #create by variable
	for(i in 1:length(fns)){ #loop thru by groups
		cat(t<-paste(strtrim(paste('--',fns[i],paste(rep('-',80),collapse='')),80),'\n')) #create line break
		q<-list(n[,2]) #create a list of group names
		names(q)<-names(n)[2] #assign by name to q
		z1<-aggregate(y,q,fns[i],na.rm=TRUE) #calculate stats
		zz1<-round(z1[,2:ncol(z1)],3) #round stats to 3 decimal places
		g<-z1[,1,,drop=FALSE] #select group variable
		z1<-cbind(g,zz1) #bind group variable with selected variables
		z2<-data.frame(t(apply(z1[,-1],1,function(x){ #calculate stats on stats
			z2<-c(nobs(x),min(x,na.rm=TRUE),max(x,na.rm=TRUE),
			mean(x,na.rm=TRUE),median(x,na.rm=TRUE),sd(x,na.rm=TRUE),cv(x))
		    names(z2)<-c('nobs','min','max','mean',
		    'median','sd','cv') #create row names
			z2<-round(z2,3) #round elements to 3 decimal places
			})))
		z<-cbind(z1,z2) #bind col stats with summary stats
		print(z) #print to console
		} #end loop thru groups
	} #end summary table w/ groups
return(z)
} #end function

"tau" <-
function(y,prior){

z<-matrix(0,nrow(y),ncol(y)) #create blank matrix
N<-sum(y)
ccr<-sum(diag(y))
n<-apply(y,1,sum)

num<-ccr-sum(prior*n)
den<-N-sum(prior*n)
tau<-num/den
tau[tau<0]<-0
return(tau)		
}

"univar.outliers" <-
function(x,var='',by='',outfile='',sd.limit=3){

if(!var==''){
	y<-subset(x,select=eval(parse(text=var))) #select variables to summarize
	}
else{y<-x}

if(by==''){ #sd outliers w/o groups
	t1<-scale(y) #calculate sd's
	t2<-t1>=sd.limit | t1<=-sd.limit #identify extreme values
	t3<-t1[apply(t2,1,FUN='any',na.rm=TRUE),,drop=FALSE] #select rows with extremes
	t4<-t3>=sd.limit | t3<=-sd.limit #identify extreme values
	t4<-t3[,apply(t4,2,FUN='any',na.rm=TRUE),drop=FALSE] #select cols with extremes
	t4[t4<sd.limit & t4>-sd.limit]<-NA #replace all other values with NA
	z<-round(t4,3) #round to 3 decimal places
	if(!outfile==''){ #write table to outfile
		write.table(z,file=paste(outfile,'.csv',sep=''),quote=FALSE,sep=',')
		} #end save outfile
	} #end sd outliers w/o groups

else{ #sd outliers w/ groups
	write('',file=outfile) #empty outfile if it exists
	n<-by.names(x,by) #create by variable
	z2<-cbind(n,y) #bind with selected variables
	groups<-levels(z2[,2]) #create vector of group names	
	for(j in 1:length(groups)){ #loop thru groups
		cat(t<-paste(strtrim(paste('--',groups[j],paste(rep('-',80),collapse='')),80),'\n')) #create line break
		t1<-z2[z2[,2]==groups[j],]  #select records for group j
		t2<-subset(t1,select=eval(parse(text=var))) #select variables
		t3<-scale(t2) #calculate sd's
		t4<-t3>=sd.limit | t3<=-sd.limit #identify extreme values
 		t5<-t3[apply(t4,1,FUN='any',na.rm=TRUE),,drop=FALSE] #select rows with extremes
		t6<-t5>=sd.limit | t5<=-sd.limit #identify extreme values
 		t6<-t5[,apply(t6,2,FUN='any',na.rm=TRUE),drop=FALSE] #select cols with extremes
		t6[t6<sd.limit & t6>-sd.limit]<-NA #replace all other values with NA
		z<-round(t6,3) #round to 3 decimal places
		print(z)
		if(!outfile==''){ #write table to outfile
			write.table(z,file=paste(outfile,'.csv',sep=''),quote=FALSE,sep=',')
			} #end save outfile
		} #end loop thru groups
	} #end sd outliers w/groups
return(z)
} #end function	

