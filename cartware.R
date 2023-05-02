
## Open Source written by Kevin McGarigal of University of Massachuestts Amerst

"cart.init" <-
function()
{
# Initialize windows for carts

	if(length(dev.list()) < 2) two.windows()
	dev.set(2)
	par(mar=rep(4,4))			# Set reasonable margins for tree

}

"two.windows" <-
function(d = '')
{
# Make two nice windows at right of screen
# Positioning depends on device name

	if(!exists('device')) device <- 'clemmys'
	device = paste(device,d,sep='')

# Default: my laptop (1400x1050), gdi interface
	if(device == 'clemmys') {
		windows(5.8,4.55,ypos=484)
		windows(5.8,4.55)
		}
		
# my laptop, sdi interface (used by turtledisp workspace)
	else if (device == 'clemmys.sdi') {
		windows(5.9,4.7,xpos=818,ypos=515)
		windows(5.9,4.7,xpos=818)
		}
	
# workstation (1280x1024)
	else if (device == 'goshawk') {
		windows(5.25,4.4,xpos=755,ypos=470)
		windows(5.25,4.4,xpos=755)
		}
	
	else if (device == 'goshawk.sdi') {
		windows(5.25,4.55,xpos=755,ypos=500)
		windows(5.25,4.55,xpos=755)
		}
# default device (small screen; e.g. cluster machines)		
	else {
		windows(4.3,3.1,ypos=344)
		windows(4.3,3.1)
		}		
}

"cart" <-
function(formula, data, ..., control=usual, pick = TRUE, prune = 0, textsize = 0.75, smooth = 1, margin=.06)
{
#Convenience shell for rpart function
#B. Compton, 2001; R version 13 May 2004
#Modified for root-only trees, 15 Dec 2004
#Opens its own windows, 25 Jul 2005
#Returns pruned tree if pick=TRUE, 2 Dec 2005

	if(length(dev.list()) < 2) cart.init()	# open windows if not yet open
	z <- rpart(formula, data, control=control, ...)
	if(smooth > 1) {
		q <- z$cptable
		for(i in 2:smooth) {
			q <- q + rpart(formula, data, control=control, ...)$cptable
		}
		z$cptable <- q/smooth
	}
	if(prune != 0)
		z <- prune.rpart(z, prune)
	if (t <- dim(z$cptable)[1] > 1) plotcp2(z)
	if (pick & t)
		z <- prune(z,pick.tree(z))
	rplot(z, textsize, margin=margin)	#	plot(z)
	z$data <- substitute(data)	# Need this for factors in cart2arc
	z
}

"text.rpart2" <-
function(x, splits = TRUE, label = "yval", FUN = text, all = FALSE, pretty = NULL, 
digits = .Options$digits - 3, use.b = TRUE, use.n = FALSE, fancy = FALSE, fwidth = 0.8, 
fheight = 0.8, ...)
{
#Modified from text.rpart to match De'ath and Fabricius 2000 output for class
#by Bradley W. Compton, 6 Dec 2000
# 30 Nov 2004 to use new rates.rpart correctly

if(!inherits(x, "rpart")) stop("Not legitimate rpart")
frame <- x$frame
col <- names(frame)
method <- x$method
ylevels <- attr(x, "ylevels")
if(!is.null(ylevels <- attr(x, "ylevels")))
col <- c(col, ylevels)
if(is.na(match(label, col)))
stop("Label must be a column label of the frame component of the tree"
)
cxy <- par("cxy")#character width and height
if(!is.null(srt <- list(...)$srt) && srt == 90)
cxy <- rev(cxy)
xy <- rpartco(x)
node <- as.numeric(row.names(x$frame))
is.left <- (node %% 2 == 0)#left hand sons
node.left <- node[is.left]
parent <- match(node.left/2, node)#Put left splits at the parent node
if(splits) {
left.child <- match(2 * node, node)
right.child <- match(node * 2 + 1, node)
rows <- labels(x, pretty = pretty)
if(fancy) {
## put split labels on branches instead of nodes
xytmp <- rpart.branch(x = xy$x, y = xy$y, node = node)
leftptx <- (xytmp$x[2,  ] + xytmp$x[1,  ])/2
leftpty <- (xytmp$y[2,  ] + xytmp$y[1,  ])/2
rightptx <- (xytmp$x[3,  ] + xytmp$x[4,  ])/2
rightpty <- (xytmp$y[3,  ] + xytmp$y[4,  ])/2
FUN(leftptx, leftpty + 0.52 * cxy[2], rows[left.child[!is.na(left.child)]], ...)
FUN(rightptx, rightpty - 0.52 * cxy[2], rows[right.child[!is.na(right.child)]], ...)
}
else FUN(xy$x, xy$y + 0.5 * cxy[2], rows[left.child], ...)
}
leaves <- if(all) rep(T, nrow(frame)) else frame$var == "<leaf>"
if(method == "class") {
if(label == "yval")
stat <- ylevels[frame$yval[leaves]]
else if(!is.na(lev <- match(label, ylevels)))
stat <- format(signif(frame$yprob[leaves, lev], digits
 = digits))
else if(label == "yprob") {
sub <- matrix(c(1:length(frame$yval), frame$yval), nrow
 = length(frame$yval))
stat <- format(signif(frame$yprob[sub][leaves], digits
 = digits))
}
else stat <- format(signif(frame[leaves, label], digits = 
digits))

r <- rates.rpart(x)
if(use.b) stat <- paste("\n", stat, "\n", (round(r[,1]/rowSums(r), 2)), "\n", "(", rowSums(r), ")", sep = "")
#if(use.b) stat <- paste("\n", stat, "\n", (round(apply(r, 1, max)/rowSums(r), 2)), "\n", "(", rowSums(r), ")", sep = "")
#stat <- paste(stat, "\n", "(", rowSums(r), ")", sep = "")
if(use.n) stat <- paste(stat,"\n", apply(r, 1, paste, collapse = "/"), sep = "")
}

else if(method == "anova") {
stat <- format(signif(frame[leaves, label], digits = digits))
	if(use.n)
		stat <- paste(stat, "\n", "(", frame$n[leaves], ")", sep = "")
	}
	else if(method == "poisson" | method == "exp") {
		stat <- format(signif(frame[leaves, label], digits = digits))
	if(use.n) {
		stat <- paste(stat, "\n", frame$yval2[leaves], "/", 
		frame$n[leaves], sep = "")
	}
}

oval <- function(middlex, middley, a, b)
{
theta <- seq(0, 2 * pi, pi/30)
newx <- middlex + a * cos(theta)
newy <- middley + b * sin(theta)
polygon(newx, newy, border = T, col = 0)
polygon(newx, newy, border = T, density = 0)
}
rectangle <- function(middlex, middley, a, b)
{
newx <- middlex + c(a, a,  - a,  - a)
newy <- middley + c(b,  - b,  - b, b)
polygon(newx, newy, border = T, col = 0)
polygon(newx, newy, border = T, density = 0)
}
if(fancy) {
## find maximum length of stat
maxlen <- max(string.bounding.box(stat)$columns) + 1
maxht <- max(string.bounding.box(stat)$rows) + 1
if(fwidth < 1)
a.length <- fwidth * cxy[1] * maxlen
else a.length <- fwidth * cxy[1]
if(fheight < 1)
b.length <- fheight * cxy[2] * maxht
else b.length <- fheight * cxy[2]
for(i in parent)
oval(xy$x[i], xy$y[i], a = (sqrt(2) * a.length)/2, b = (
sqrt(2) * b.length)/2)
child <- match(node[frame$var == "<leaf>"], node)
for(i in child)
rectangle(xy$x[i], xy$y[i], a = a.length/2, b = 
b.length/2)
}
#if FUN=text then adj=1 puts the split label to the left of the
#    split rather than centered
#Allow labels at all or just leaf nodes
## stick values on nodes

if(fancy)
	FUN(xy$x[leaves], xy$y[leaves] + 0.5 * cxy[2], stat, ...)
else FUN(xy$x[leaves], xy$y[leaves] - 0.5 * cxy[2], stat, adj = 0.5, 
...)
invisible()
}

"flop.rpart" <-
function(x, n = 1, m = 1)
{
# flop.rpart(rpart object)
# Flop branches of rpart object x such that splits are always > (or :), never <
# B. Compton, 2 March 2003

#return(x)

	if(n == 1) {
		y <- x
		x <- x$frame
	}
	i <- x[as.numeric(row.names(x)) == n,  ]
	row.names(i) <- as.character(m)
	if(i[1] == "<leaf>")
		i
	else {
		q <- c(2 * n, 2 * n + 1)	# Existing node #s
		r <- c(2 * m, 2 * m + 1)	# New node #s
print(i$splits)
		if("<" == substring(i$splits[1], 1, 1)) {
			i$splits[1] <- paste(">", substring(i$splits[1], 2), 
				sep = "")
			i$splits[2] <- paste("<", substring(i$splits[2], 2), 
				sep = "")
			q <- rev(q)
		}
		z <- rbind(flop.rpart(x, q[1], r[1]), flop.rpart(x, q[2], r[2])
			)
		z <- rbind(i, z)	# if within recursion, return new node
		if(n != 1)
			z
		else {
			y$frame <- z
			y
		}
	}
}

"pick.tree" <-
function(a) 
{
# pick.tree
# Pick best tree from rpart object using 1 SE rule
# Returns number of leaves in global leaves
# B. Compton, 12 May 2004, 21 May 2004

leaves <<- 1+as.numeric(a$cptable[,2][i <- a$cptable[,4]<=min(a$cptable[,4]+a$cptable[,5])][1])
a$cptable[,1][i][1] + .000001   # Add a little for rounding errors

}

"splits.rpart" <-
function(z) {
# splits.rpart
# B. Compton, 19 Apr 2004
# Gives results from old z$frame$splits

z <- labels(z, collapse = FALSE)
z[a <- apply(z, 2, "==", "<leaf>")] <- ""
a <- (!a) & apply(z, 2, substring, 1, 1) != "<" & apply(z, 2, substring, 1, 1) != ">"
z[a] <- paste(":", z[a], sep="")
z[!a] <- paste( apply(z, 2, substr, 1, 1),apply(z, 2, substring, 3, last=100), sep="")[!a]
colnames(z) <- c("cutleft", "cutright")
z 
}

"cart2arc" <-
function(a, filename = "zzz.aml", result = "z", q, i = 1, l = 0, avoid = "", mask="", source = "")
{
# cart2arc: Create Arc/Info Grid AML to apply results of rpart (CART) analysis in R
# Parameters:
# a			- an rpart object
# filename	- filename to route output to
# result	- name of Arc grid to create as result
# source	- path to source grids
# avoid		- include avoid = clause to restrict to certain areas, e.g., avoid = "ocean == 1" 
#			  will give zeros where grid ocean is 1
# mask		- include mask = grid to set mask (in Arc), giving NODATA for masked out cells
# Global class.code contains id -> gridcode substitutions
#	class.code is an n x 1 character or numeric matrix; row names are factors in the Y
#   variable, and values are corresponding arc grid codes
#	Set it like this:
#		class.code <- matrix(c(0,1),2,1)
#		row.names(class.code) <- c("absent","present")
#	set class.code <- as.null() if you don't want to use it.
# Global var.file contains variable -> arc grid file substitutions
#	var.file is an n x 1 character matrix; row names are variables, values are filenames
#	set var.file <- as.null() if you don't want to use it.
# N.B. When mapping presence/absence, 'present' must be the first factor in your data, because
# if 'absent' is first, present will be pushed to the ELSE clause, and cover types for which 
# there is no data will be mapped as present.
# Maximum of 26 levels of factor variables
# Bradley W. Compton
# 21 Dec 2000
# Modified to work in R, 19 Apr 2004
# Modified more to work with stupid Arc/Grid, 13 May 2004
# Modified for easier class.code and var.file, new mask & avoid, etc. 19 Feb 2005

 	if(dim(a$cptable)[1] == 1)
		return(warning("   cart2arc fed empty tree; no results produced"))

	old <- options(width = 255)
	on.exit(options(old))
	if(l == 0) {											# if top level, prepare data
		q <- data.frame(as.numeric(row.names(a$frame)), gsub("[\n]","_",I(as.character(a$frame$var))),
			splits.rpart(a), attr(a, "ylevels")[a$frame$yval])
		dimnames(q)[2] <- list(c("branch", "var", "cutleft", "cutright", "y"))
		cat("/*", filename, file = filename, fill = T, append = F)
		cat("/* Automatically created by cart2arc,", date(), file = filename, "\n", fill = TRUE, append = TRUE)
		if(0!=nchar(mask)) 
			cat("setmask",mask,"\n")
		cat("docell", file = filename, fill = TRUE, append = TRUE)
		if(0!=nchar(avoid))
			cat("if (", avoid, ")\n   ", result, " = 0\nelse", file = filename, 
			fill = TRUE, append = TRUE, sep = "")
	}
	j <- (1:(dim(q)[1]))[q$branch == i]						#index into q
	if(q$var[j] == "<leaf>") {								# if terminal leaf
		t <- paste(rep(" ", (l + 1) * 3), collapse = "")	#    write out assignment
		k <- as.character(q$y[j])
		if (!is.null(class.code)) k <- class.code[k,]
		cat(t, result, " = ", k, "\t/* Class = ", as.character(q$y[j]), file = 	filename, fill = TRUE, append = TRUE, sep = "")
	}
	else {													# else, write if-statement for node
		t <- paste(rep(" ", (l + 1) * 3), collapse = "")	#  tabs
		s <- as.character(q$cutleft[j])
		e <- ""
		if(substring(s, 1, 1) == ":") {						#  if variable is factor
			s <- attr(eval(a$data)[, as.character(q$var[j])], "levels")[match(
				substring(s, 2:nchar(s), 2:nchar(s)), letters)]
			s <- paste(", 'value in {", paste(as.numeric(s), collapse = ", "
				), "}')", sep = "")
			e <- "test("
		}
		else s <- paste(substring(s, 1, 1), substring(s, 2))

		if(is.null(var.file) || is.na(v <- as.character(var.file[as.character(q$var[j]),])))
			v <- paste(source, as.character(q$var[j]), collapse = "", sep = "")
		cat(t, "if (", e, v, " ", s, ")", file = filename, fill = TRUE, append = TRUE, sep = "")
		cart2arc(a, filename, result, q, q$branch[j] * 2, l + 1, avoid=avoid, source=source)		#  recurse for left split
		cat(t, "else", file = filename, fill = TRUE, append = TRUE, sep = "")
		cart2arc(a, filename, result, q, q$branch[j] * 2 + 1, l + 1, avoid=avoid, source=source)	#  recurse for right split
	}
	if(l == 0)
		cat("end", file = filename, fill = TRUE, append = TRUE)
		if(0!=nchar(mask)) 
			cat("setmask off\n")
}

"confuse" <-
function(a)
{
#Return confusion matrix for rpart object
#Completely rewritten version
#B. Compton, 11 Nov 2004
#I think this will work correctly, including with priors and a loss matrix
#Transposed to match SAS - rows are observed, columns are predicted

	y <- attr(a, 'ylevels')
	x <- cbind(a$y, a$frame$yval[a$where])
	z <- matrix(0,length(y),length(y))
	for(i in 1:dim(x)[1]) z[x[i,1],x[i,2]] <- z[x[i,1],x[i,2]] + 1
	rownames(z) <- y
	colnames(z) <- y
	z
}

"rplot" <-
function(z, textsize = 0.75, misclass = FALSE, digits = 5, quiet = F, margin=.06, use.b=TRUE, use.n=FALSE, ...)
{
# Cover function to plot rpart objects
# B. Compton; last update 2 Mar 2003


#	z <- flop.rpart(z)	# Make trees always use > splits
#	r <- rates.rpart(z)

	if(z$method == "anova") use.n=TRUE
	plot(z, margin = margin, ...)
	text.rpart2(z, use.b=use.b, use.n=use.n, cex = textsize, digits = digits)
	if(!quiet) {
#	legend(1, 1, paste(names(colSums(r)), colSums(r), sep = ": "), bty = "n")
		if(z$method == "anova") {
			cat("Tree has", rev(z$cptable[, 2])[1] + 1, "leaves.  ")
			cat("R-squared =", round((1 - rev(z$cptable[, 3]))[1], 
				digits = 4), "\n")
		}
		else {
			c <- confuse(z)
			m <- round(c(100*max(rowSums(c))/sum(c), 100*sum(diag(c))/sum(c), sum(diag(c)), sum(c)))
			q <- paste("Correct classification rate: Null = ", m[1], "%, Model = ", m[2], "% (", m[3], "/", m[4], ")", sep = "")
			
			if(misclass)
				legend(locator(1), q, bty = "n")
				
			cat(c(q, "\n"))
			cat(c(kappatau(c, z$parms$prior), "\n"))

			# Get number of leaves.  An error in prune allows a CP to be repeated sometimes, so we need the first close match.
			t <- z$cptable[dim(z$cptable)[1],1]
 			q <- vector(length = dim(z$cptable)[1])
			for(i in 1:dim(z$cptable)[1]) q[i] <- identical(all.equal(z$cptable[i,1], t),TRUE)
			l <- z$cptable[q,2][1]

			cat("Leaves = ",l + 1,"\n")
			cat("\n")
			cat("Confusion matrix (rows = observed, cols = predicted)\n")
			print(c)
		}
	}
}

"monte.cart" <-
function(model, data, n = 1000, size, ifplot=FALSE, ifcartplot=TRUE, control=usual, parms,...)
{
# monte.cart
# Monte Carlo CART resampling.  Gives P value for rpart tree.
# What proportion of r^2's from trees based on randomly permuted data
# are as extreme as the r^2 from my tree?
# This version throws out any trees that don't match the size of my tree.
# Result[1] = P-value;
# Result[2:n] = vector of r/CCRs
# B. Compton, 14-20 Feb 2003
# 28 May 2004: Modified to work for classification trees too, using ccr.
# 31 May 2004: Only permutes variables in model, not whole dataset
# 6-7 Jun 2004: Fixed bug in calculating CCR
# 11 Nov 2004: Permute first response variable only to preserve correlation structure
# 11 Nov 2004: Give P < x to incorporate number of reps in statistic
# 11 Nov 2004: Give error if size unavailable; set minsplit & minbucket!
# 23 Nov 2004: Use same default (usual) for control as cart; allow specifying control
# 4 Feb 2005: Select subset of variables with all.vars
# 14 Feb 2005: Gives warning if can't find many trees of proper size
# 23 Mar 2005: Allow passing parms (for loss matrix), add fudge factor to get correct tree size, 
# 	and calculate CCR directly, as opposed to improvement over null.
# 30 Mar 2005: Fix silly bug in permuted CCRs


	fudge <- 0.0001 	# Fudge factor for complexity
	d <- data[p <- all.vars(model, max.names=10000)]
	cat("\n\n",title <- paste("monte.cart for ",p[1],", size = ",size, sep=""),"\n")
	cat("   Model = ")
	print(model)
	z <- rep(0, n)
	control$cp <- 1e-005
	control$maxsurrogate <- 0
	control$xval <- 0
	control$maxdepth <- size

	if(missing(parms)) 
		f <- rpart(model, data = d, control = control, ...)
	else 
		f <- rpart(model, data = d, control = control, parms = parms, ...)

	if (f$method == "class") {
		r <- f$cp[f$cp[, 2] == (size - 1), 1]
		if(length(r) != 0) {
	 		r <- confuse(prune(f,r + fudge))
	 		r <- sum(diag(r))/sum(r)			# Correct classification rate
			t <- "CCR"
			}
		}
	else {
		t <- "R-squared"
		r <- 1 - (f$cp[f$cp[, 2] == (size - 1), 3])
		}
		if (length(r) == 0)
			stop(paste("No tree of size",size,"available"))
		cat(paste("  ",t,"=", round(r, digits = 4)), "\n")

	if(ifcartplot)
		rplot(prune(f,f$cp[f$cp[, 2] == (size - 1), 1] + fudge), textsize=0.75)
		# Add fudge factor to complexity to make sure we don't get too many leaves.  Ugh.

	i <- 1
	j <- 0
	while(i <= n) {
		d[,1] <- sample(d[,1]) 				# Permute first variable
		if(missing(parms))
			f <- rpart(model, data = d, control = control, ...)
		else	
			f <- rpart(model, data = d, control = control, parms = parms,...)

		if (f$method == "class") {
			q <- f$cp[f$cp[, 2] == (size - 1), 1] 
			if(length(q) != 0) {
				q <- confuse(prune(f, q + fudge))
 				z[i] <- sum(diag(q))/sum(q)			# Correct classification rate
 				}
 			else z[i] <- NA
 			}
		else
			z[i] <- 1 - (rev(f$cp[f$cp[, 2] == (size - 1), 3]))[1]
		if(!is.na(z[i])) 
			i <- i + 1
		else
			j <- j + 1
		if((j < 1e6) & (j/i >= 100)) {
			cat(paste("Warning: Having trouble finding ",size,"-leaf trees (",i-1," successes, ",j," attempts).\n",sep = ""))
			cat("   Press Esc to give up.\n")
			j <- 1e6
		}
	}
	cat("\nPermutation test results:\n")
	P <- round((sum(z >= r) + 1)/n, digits = 4)
	if (P > 1) 
		cat("   P = 1\n")
	else
		cat(paste("   P < ", P, sep = ""), "\n")
	d <- density(z, bw = bw.SJ(z, method = "dpi"))
	cat(paste("   kernel-based P = ",round(sum(d$y[d$x >= r])/sum(d$y), digits=4), sep = ""), "\n")

	if (ifplot) {
		title <- paste(title," (P < ",P,")",sep="")
		dev.set(3)	
		plot(d,xlim=c(min(d$x,r),max(d$x,r)),main=title,xlab=paste("n =",n))
		points(r,0,pch=24,bg=1,cex=1.5)
		dev.set(2)
	}
	c(P,z)
}

"loss" <-
function(cost)
{
# Create 2x2 loss matrix from cost

if (length(cost) == 1)
		cost <- c(cost,1)
list(loss = t(matrix(c(0,cost,0),2,2)))
}

"plotcp2" <-
function(...)
{
	dev.set(3)
	plotcp(...)
	dev.set(2)
	return()
}

"rates.rpart" <-
function(a) 
{
# Get misclassification rates for rpart object
# B. Compton, 13 May 2004
# Rewritten, 16 Nov 2004
# Modified for root trees, 15 Dec 2004

x <- a$frame$yval2[a$frame$var == "<leaf>",]
if(is.matrix(x)) {
	z <- matrix(0,dim(x)[1], 2)
	for(i in 1:dim(x)[1]) 
		z[i,1] <- x[i,x[i,1]+1]
	z[,2] <- rowSums(x[,2:((dim(x)[2]-1)/2+1)]) - z[,1]
	}
else z <- matrix(x[c(3,2)],1,2)
colnames(z) <- c('right', 'wrong')
z

}

"kappatau" <-
function(c, p)
{
# kappatau
# Compute Kappa or Tau statistic for confusion matrix.  Returns as global kt.
# Reference: McGarigal, K, S. Cushman, and S. Stafford. 2000. Multivariate statistics
# 	 for wildlife and ecology research. Springer-Verlag, New York. Pages 165-166.
# B. Compton, 15 Nov 2004, 12 Aug 2005, 29 Nov 2005, 8 May 2006

	b <- 0 != rowSums(c)+colSums(c)		# remove any unused classes
 	c <- c[b,b]
	p <- p[b]
	if (all(p == (rowSums(c)/sum(c)))) {
		# Priors are proportional to group sizes (default), so do kappa
		e <- sum(rowSums(c)*colSums(c))/sum(c)
		kt <<- (sum(diag(c)) - e) / (sum(c) - e)
	
		cat ('Kappa = ',round(kt,3), '\n')
		}
	else {
		# Priors specified, so do tau
		e <- sum(p * rowSums(c))
		kt <<- (sum(diag(c)) - e) / (sum(c) - e)
	
		cat ('Tau = ',round(kt,3), '\n')
		}
}

"usual" <-
structure(list(minsplit = 8, minbucket = 5, cp = 0, maxcompete = 100, 
    maxsurrogate = 100, usesurrogate = 2, surrogatestyle = 0, 
    xval = 10, maxdepth = 30), .Names = c("minsplit", "minbucket", 
"cp", "maxcompete", "maxsurrogate", "usesurrogate", "surrogatestyle", 
"xval", "maxdepth"))
"rpartco" <-
function(tree, parms =  paste(".rpart.parms", dev.cur(), sep = "."))
    {

    frame <- tree$frame
    node <- as.numeric(row.names(frame))
    depth <- tree.depth(node)
    is.leaf <- (frame$var == '<leaf>')
    if (exists(parms, envir=.GlobalEnv)) {
parms <- get(parms, envir=.GlobalEnv)
uniform <- parms$uniform
nspace <-parms$nspace
minbranch <- parms$minbranch
}
    else {
uniform <- FALSE
nspace <- -1
minbranch <- .3
        }

    if(uniform) y <- (1 + max(depth) -depth) / max(depth,4)
    else {                    #make y- (parent y) = change in deviance
y <- dev <- frame$dev
        temp <- split(seq(node), depth)     #depth 0 nodes, then 1, then ...
        parent <- match(floor(node/2), node)
        sibling <- match(ifelse(node %% 2, node - 1, node + 1), node)

# assign the depths
        for(i in temp[-1]) {
    temp2 <- dev[parent[i]] - (dev[i] + dev[sibling[i]])
            y[i] <- y[parent[i]] - temp2
    }
#
# For some problems, classification & loss matrices in particular
#   the gain from a split may be 0.  This is ugly on the plot.
# Hence the "fudge" factor of  .3* the average step
#
fudge <-  minbranch * diff(range(y)) / max(depth)
        for(i in temp[-1]) {
    temp2 <- dev[parent[i]] - (dev[i] + dev[sibling[i]])
    haskids <- !(is.leaf[i] & is.leaf[sibling[i]])
    y[i] <- y[parent[i]] - ifelse(temp2<=fudge & haskids, fudge, temp2)
    }
y <- y / (max(y))
        }

    # Now compute the x coordinates, by spacing out the leaves and then
    #   filling in
    x   <-  double(length(node))         #allocate, then fill it in below
    x[is.leaf] <- seq(sum(is.leaf))      # leaves at 1, 2, 3, ....
    left.child <- match(node * 2, node)
    right.child <- match(node * 2 + 1, node)

    # temp is a list of non-is.leaf, by depth
    temp <- split(seq(node)[!is.leaf], depth[!is.leaf])
    for(i in rev(temp))
            x[i] <- 0.5 * (x[left.child[i]] + x[right.child[i]])

    if (nspace < 0) return(list(x=x, y=y))

    #
    # Now we get fancy, and try to do overlapping
    #
    #  The basic algorithm is, at each node:
    #      1: get the left & right edges, by depth, for the left and
    #           right sons, of the x-coordinate spacing.
    #      2: find the minimal free spacing.  If this is >0, slide the
    #           right hand son over to the left
    #      3: report the left & right extents of the new tree up to the
    #           parent
    #   A way to visualize steps 1 and 2 is to imagine, for a given node,
    #      that the left son, with all its descendants, is drawn on a
    #      slab of wood.  The left & right edges, per level, give the
    #      width of this board.  (The board is not a rectangle, it has
    #      'stair step' edges). Do the same for the right son.  Now
    #      insert some spacers, one per level, and slide right hand
    #      board over until they touch.  Glue the boards and spacer
    #      together at that point.
    #
    #  If a node has children, its 'space' is considered to extend left
    #    and right by the amount "nspace", which accounts for space
    #    used by the arcs from this node to its children.  For
    #    horseshoe connections nspace usually is 1.
    #
    #  To make it global for a recursive function, the x coordinate list
    #    is written into frame 0.
    #
    compress <- function(me, depth) {
        lson <- me +1
x <- x
if (is.leaf[lson]) left <- list(left=x[lson], right=x[lson],
depth=depth+1, sons=lson)
        else               left <- compress(me+1, depth+1)

        rson <- me + 1 + length(left$sons)        #index of right son
if (is.leaf[rson]) right<- list(left=x[rson], right=x[rson],
depth=depth+1, sons=rson)
else               right<- compress(rson, depth+1)

maxd <- max(left$depth, right$depth) - depth
        mind <- min(left$depth, right$depth) - depth

# Find the smallest distance between the two subtrees
#   But only over depths that they have in common
# 1 is a minimum distance allowed
slide <- min(right$left[1:mind] - left$right[1:mind]) -1
if (slide >0) { # slide the right hand node to the left
    x[right$sons] <- x[right$sons] - slide;
    x[me] <- (x[right$sons[1]] + x[left$sons[1]])/2
#    assign("x", x)
            x <<- x
    }
else slide <- 0

# report back
        if (left$depth > right$depth) {
    templ <- left$left
            tempr <- left$right
            tempr[1:mind] <- pmax(tempr[1:mind], right$right -slide)
    }
        else {
    templ <- right$left  - slide
    tempr <- right$right - slide
    templ[1:mind] <- pmin(templ[1:mind], left$left)
    }

list(left = c(x[me]- nspace*(x[me] -x[lson]), templ),
     right= c(x[me]- nspace*(x[me] -x[rson]), tempr),
     depth= maxd+ depth, sons=c(me, left$sons, right$sons))
}
#    assign('compress', compress)
#    assign('x', x)
#    assign('is.leaf', is.leaf)
#    assign('nspace', nspace)

#    temp <-
    compress(1, 1)
#    x <- get('x')
#    remove(c('compress', 'x', 'is.leaf', 'nspace'))
    list(x = x, y = y)
}

"tree.depth" <-
function (nodes)
{
    depth <- floor(log(nodes, base = 2) + 1e-7)
    as.vector(depth - min(depth))
}

"var.importance" <-
function(x, chatter = FALSE) {
# importance
# Give variable importances for classification or regression tree x
# Note that less-important variables will be omitted if 
# rpart.control(maxsurrogate = ...) is small.
# Source: Breiman et al. (1993), page 146-150.
# B. Compton, 17-18 and 22-23 Nov 2005
# to do: matches CART for regression trees; close but not perfect for classification trees
#	still doesn't deal with priors and costs
# Subfunctions: improve, ssd, gini


	f <- x$frame									# Frame for rpart object
	g <- data.frame(0,x$splits[,c('improve','index','ncat')],row.names=NULL)	# All primary, alternate, and surrogate splitters
	g[,1] <- row.names(x$splits)
	names(g)[1] <- 'variable'
	m <- g[FALSE,c('variable','improve')]			# Result to build
	for(i in 1:dim(f)[1]) {							# For each row of frame,
		if(f[i,'var'] != '<leaf>') {				#	Only interested in nodes
			n <- as.numeric(row.names(f)[i])		#	Node number
			k <- sum(f[(1:i)[-i],'ncompete']) + sum(f[(1:i)[-i],'nsurrogate']) + sum(f[(1:i)[-i],'var'] != '<leaf>') + 1
			p <- g[k,c('variable','improve')]
#	print(g[k,])
#	print(n)
				t <- (snip.rpart(x,n))					#	Tree snipped to this node

				b <- t$where == (1:dim(t$frame)[1])[row.names(t$frame) == n]	#	Rows in this node
#	print(table(b))
				v <- eval(x$data)[,g[k,'variable']]				#		Values of surrogate variable
				l <- x$y[b & v < g[k,'index']]				#			Left child
				r <- x$y[b & v >= g[k,'index']]				#			Right child
#	catn('split=',g[k,'index'])
#	catn('v=',v)
#	catn('l=',l)
#	catn('r=',r)
				p2 <- improve(x$y[b], l, r, x$method, x$parms$loss, as.vector(x$parms$prior), x$y)
#return()
			#	p2 <- p2 / x$frame[i,'n']			#	Adjust by n to match improvement in CART

			if(x$method == 'anova') p[2] <- p[2] * x$frame[i,'dev']			# 	Adjust by deviance for regression trees
# To match CART for regression trees, following line needs to divide by total n; but
# it seems that it should use row's n for classification trees or something....
#			p[2] <- p[2] / x$frame[1,'n']			#	Adjust by n to match improvement in CART [matches CART for regression trees]
			p[2] <- p[2] / x$frame[i,'n']			#	Adjust by n to match improvement in CART [correct for classification trees]
			
			m <- rbind(m, p)						#	Get improvement of primary splitter
if(chatter) cat(paste('node ',i,', primary ',p[1],', improve = ',round(p[2],8),'\n', sep=''))
if(chatter) cat(paste('node ',i,', primary ',p[1],', calc. i = ',round(p2,8),'\n', sep=''))

			# Now get improvement for each surrogate
			t <- (snip.rpart(x,n))					#	Tree snipped to this node
			b <- t$where == (1:dim(t$frame)[1])[row.names(t$frame) == n]	#	Rows in this node
			s <- f[i,'nsurrogate']					#	Number of surrogates at this node
			s <- sum(f[1:i,'ncompete']) + sum(f[(1:i)[-i],'nsurrogate']) + (1:(s+1))[-(s+1)] + sum(f[1:i,'var'] != '<leaf>')
			for(j in s) {										#	For each surrogate at this node (j points to x$splits),
				v <- eval(x$data)[,g[j,'variable']]				#		Values of surrogate variable
				if(abs(g[j,'ncat']) > 1) {						#		If categorical predictor,
					c <- x$csplit[g[j,'index'],]				#			Categories
					l <- x$y[b & v %in% (levels(v))[c == 1]]	#			Left child
					r <- x$y[b & v %in% (levels(v))[c == 3]]	#			Right child
				}
				else {											#		else, continuous predictor
					l <- x$y[b & v < g[j,'index']]				#			Left child
					r <- x$y[b & v >= g[j,'index']]				#			Right child
				}
				z <- improve(x$y[b], l, r, x$method, x$parms$loss, as.vector(x$parms$prior), x$y)
				if(x$method == 'anova') z <- z / x$frame[1,'n']	#		If regression tree, adjust deviance by sample size
				m <- rbind(m,list(g[j,'variable'],z))
### if(chatter) cat(paste('node ',i,', surrogate [',j,'] ',g[j,'variable'],', improve = ',round(z,8),'\n', sep=''))
			}	
		}
	}
	m <<- m
	z <- aggregate(m$improve, by=list(m$variable), FUN='sum')	# Sum delta-I's across variables
	names(z) <- c('variable','importance')
	z$importance <- 100 * z$importance / max(z$importance)		# Relative scaling
	z$importance <- round(z$importance, 2)
	z <- z[order(z$importance, decreasing = TRUE),]
	row.names(z) <- 1:dim(z)[1]
	z
}

"improve" <-
function(n,l,r,t,c,p, y) {
# Calculate delta-I for values in node, left child, and right child
# Use deviance for regression trees (t = 'anova') and gini for classification trees
# c = costs, p = priors, y = values of response variable

catn<<-function(...){cat(paste(...,'\n'))}

	n <- n[!is.na(n)]
	l <- l[!is.na(l)]
	r <- r[!is.na(r)]

	if(t == 'anova') {		# If regression tree, use deviance
		z <- ssd(n) - ssd(l) - ssd(r)
		if(is.na(z)) catn(n,l,r)
		}
	else {					# For classification trees, use Gini index

catn('length(l)=',length(l),'length(r)=',length(r),'length(n)=',length(n))

		q <- c(length(l),length(r)) / length(n)
#		z <- gini(n, c, p, y) - sum(q * c(gini(l, c, p, y), gini(r, c, p, y)))

a <- gini(n, c, p, y) 
b <- gini(l, c, p, y)
c <- gini(r, c, p, y)

z <- a - sum(q * c(b,c))
catn('a=',a,'b=',b,'c=',c,'q=',q,'improvement=',z)

	}
#nn<<-n;cc<<-c;pp<<-p;yy<<-y;ll<<-l;rr<<-r;qq<<-q
	z
}

"ssd" <-
function(x)
# Give sum of squared deviations of x
{
	sum((x - mean(x))^2)
}

"gini" <-
function(x, c = NULL, p = NULL, y)
# Gini index of x
# Use costs c and priors p, and entire response variable y
# B. Compton, 29 Nov 2005
#
# Current version (14 Dec 2005) implements priors but not costs.
# It matches Kevin's Gini/priors results, but still doesn't
# give me same improvement as rpart for primary splitters when
# priors are set (it matches when priors come from data).
# Numbers are close (within 15%), but with no clear pattern of
# disagreement.  Either Kevin & I have both misinterperted how
# the Gini index incorporates priors, or the weightings for
# the left/right splits are wrong.  Grr.


{
#    z <- rowSums(outer(unique(x),x,FUN='=='))
#    catn('gini1=',1 - sum((z / sum(z))^2))

#x<<-x;c<<-c;p<<-p;y<<-y
#catn('x=',x)
#catn('c=')
#print(c)
#catn('p=',p)
#catn('y=',y)


    g <- length(u <- unique(y)) # number of groups (get levels from response variable)
    h <- table(y) # total tree count in each group
    n <- rowSums(outer(u,x,FUN='==')) # number in each group
    s <- sum(n) # n across groups
#catn('g=',g)
#catn('u=',u)
#catn('h=',h)
#catn('n=',n)


    if(is.null(c))    c <- 1*outer(1:g,1:g,FUN='!=') # default costs = all 1
    if(is.null(p))    p <- h/sum(h) # default priors = from data
###p <- p*g  #    p <- p/min(p)  #???

#    c <- c * p    

#    print(c)                                    # express priors as costs
#g<<-g;n<<-n;s<<-s;c<<-c;p<<-p;m<<-m
#catn('gini2=',sum(c * m * t(m)))
#m <- matrix(n,g,g) / s
#1 - sum(p * n / sum(n))
#    1-sum(c * m * t(m))                         # Eq. 4.12 from Breiman et al.
#1 - sum(((n/sum(n))^2))


#catn('p=',p)
#catn('n=',n)
#catn('h=',h)

    r <- p*n/h                    # Kevin's version
    1 - sum((r/sum(r))^2)

}

"rhist" <-
function(a, digits=4) {
# rhist
# Give histogram for each leaf of classification or regression trees (rpart object a)
# B. Compton, 20 Dec 2005
# 4 Jan 2006: use silly little hats in plots


	r <- (1:dim(a$frame)[1])[a$frame[,'var'] == '<leaf>']
	v <- a$frame[r,'yval']
	
	n <- floor(length(r)^.5)
	n <- c(n, ceiling(length(r)/n))
	q <- dev.set()
	windows(2*n[2],2*n[1])
	par(mfrow=n)

	for(i in 1:length(r)) {
		if(a$method == 'class') {
			u <- sort(unique(a$y))
			b <- rowSums(outer(u, a$y[a$where == r[i]], FUN='=='))
			barplot(b, names.arg=u, col=c('gray','pink2')[(u == v[i]) + 1], xlab=bquote(hat(y) == .(v[i])))
		}
		else {
			hist(a$y[a$where == r[i]], xlim=c(min(a$y),max(a$y)), xlab=bquote(hat(y) == .(round(v[i],digits))), main=NULL, ylab=NULL)
 			lines(rep(v[i],2), c(-1,10000), col='red')
 		}
 	}
 	q <- dev.set(q)
}

