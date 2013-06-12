#!/usr/bin/env Rscript 
library("biglm")
suppressPackageStartupMessages(library("optparse"))

option_list <- list(

	make_option(c("-f", "--FileName"), action="store", type = "character", default=NULL,
				help="Fasta file to be analyzed [default = %default]"),
				
	make_option(c("-c", "--cycles"), action="store", type = "integer", default=15,
				help="Iterations before widening prediction interval [default = %default]"),

	make_option(c("-t", "--Tthreshold"), action="store", type="double", default=1e-05,
				help="Threshold value for parameter estiamte t [default = %default]"),

	make_option(c("-m", "--Mthreshold"), action="store", type="double", default=1e-03,
				help="Threshold value for parameter estiamte m [default = %default]"),

	make_option(c("-l", "--Lthreshold"), action="store", type="double", default=8e-05,
				help="Threshold value for parameter estiamte l [default = %default]"),
	
	make_option(c("-q", "--quiet"), action="store", type="character", default="F",
				help="disable all warning messages (set to T or F) [default = %default]"), 
	
	make_option(c("-o", "--output"), action="store", type="character", default=NULL,
				help="save output to a .txt file [default = %default]")
			
)
opt <- parse_args(OptionParser(option_list = option_list, 
				  description = "\n Optimization procedure that estimates evolutionary distance between two sequences (t), mean of indel length (m), and the rate of insertion and the rate of deletion per substitution (l) given raw alignment sequences in a .fasta file format"))

if (opt$quiet == "T") {
	options(warn = -1)
}

if (is.null(opt$output) == F) {
	write(c("Parameter t", "Parameter m", "Parameter l"), file = opt$output,
		  ncolumns = 3, append = FALSE, sep = "\t")
}
########################################################################################
cat("\n")
########################################################################################
start_time <- proc.time()

if (is.null(opt$FileName) == T) { 
	ExactFile<-readLines(file("example.dawg", "r"))
	system ('dawg - -o R.fasta --label=on', intern = FALSE, input = ExactFile)
} else {
	invisible(file.copy(opt$FileName, "R.fasta"))
}

#system ('perl -i -pe \'s/%r/++$i/eg\' R.fasta', intern = FALSE)
system ('./NgilaWrapper R.fasta', intern = FALSE) 

sim_values <- read.table ("sim_values.txt", col.names = c("BL_exact", "IM_exact", "L_exact", "nuisance"))

attach(sim_values)

BL <- BL_exact
M <- IM_exact
L <- L_exact

sink("Target_Values.txt")
cat(BL_exact, "\n")
cat(IM_exact, "\n")
cat(L_exact, "\n")
sink()

system('rm Output.fasta', intern = FALSE)
system('rm sim_values.txt', intern = FALSE)

closeAllConnections() 
detach(sim_values)

########################################################################################

load("predict.RData")

# Prediction Interval
w <- 0.90

mm <- matrix(data = c(BL, M, L, BL*M, BL*L, M*L, BL*M*L), ncol = 7)

t_fit <- (t_func$fit(mm))/2
t_se <- (sqrt(diag((t_func$se(mm))) + t_reg$qr$ss/t_reg$df))/2
cat("t fit is: ", t_fit, "\n")

m_fit <- m_func$fit(mm)
m_se <- sqrt(diag(m_func$se(mm)) + m_reg$qr$ss/m_reg$df)
cat("m fit is: ", m_fit, "\n")

l_fit <- l_func$fit(mm)
l_se <- sqrt(diag(l_func$se(mm)) + l_reg$qr$ss/l_reg$df)
cat("l fit is: ", l_fit, "\n")
cat("\n")

# calculate prediction interval
interval_t <- cbind(t_fit,t_fit)+cbind(t_se,-t_se)*qt((1-w)/2,t_reg$df)
interval_m <- cbind(m_fit,m_fit)+cbind(m_se,-m_se)*qt((1-w)/2,m_reg$df)
interval_l <- cbind(l_fit,l_fit)+cbind(l_se,-l_se)*qt((1-w)/2,l_reg$df)

update_t <- t_fit
update_m <- m_fit
update_l <- l_fit
update_ir <- l_fit / (2 * t_fit)

cat("interval t is: ", interval_t, "\n")
cat("interval m is: ", interval_m, "\n")
cat("interval l is: ", interval_l, "\n")
cat("\n")

########################################################################################

count_cycles  <- 1
finished <- 0
while (finished == 0) {

	T_set <- function(x) {

		data <- file("Target_Values.txt", "r")

		BL_exact <- readLines(data, n=1)
		BL_exact <- as.numeric(BL_exact)
		
		BL_experimental <- x[1]
		
		sink("Final.dawg")
		cat("Tree.Tree = (C_%r:",BL_experimental,",A_%r:",BL_experimental,");")
		cat("\n")
		cat("Root.Length = 1000")
		cat("\n")
		cat("Subst.Model = jc")
		cat("\n")
		cat("Subst.Freqs = 0.2,0.3,0.3,0.2")
		cat("\n")
		cat("Indel.Model = geo")
		cat("\n")
		cat("Indel.Params =", update_m)
		cat("\n")
		cat("Indel.Rate =", update_ir)
		cat("\n")
		cat("Sim.Reps = 100")
		cat("\n")
		sink()

		EstimateFile<-readLines(file("Final.dawg", "r"))
		system ('dawg - -o R.fasta --label=on', intern = FALSE, input = EstimateFile)
		#system ('perl -i -pe \'s/%r/++$i/eg\' R.fasta', intern = FALSE)
		system ('./NgilaWrapper R.fasta',intern = FALSE)	

		sim_values <- read.table ("sim_values.txt", 
		col.names = c("BL_final", "IM_final", "L_final", "nuisance"))
		attach(sim_values)
		
		check_t <- ((BL_exact - BL_final)^2 / BL_exact)

		system('rm Output.fasta', intern = FALSE)
		system('rm sim_values.txt', intern = FALSE)
		detach(sim_values)
		return (check_t)	
	}

	answer_t <- optimize(T_set, interval_t, tol = 0.001)
	update_t <- answer_t$minimum
	update_ir <- update_l / (2*update_t)
	check_t <- answer_t$objective

########################################################################################

	M_set <- function(x) {


		data <- file("Target_Values.txt", "r")

		IM_exact <- readLines(data, n=1)
		IM_exact <- readLines(data, n=1)
		IM_exact <- as.numeric(IM_exact)

		IM_experimental <- x[1]

		sink("Final.dawg")
		cat("Tree.Tree = (C_%r:",update_t,",A_%r:",update_t,");")
		cat("\n")
		cat("Root.Length = 1000")
		cat("\n")
		cat("Subst.Model = jc")
		cat("\n")
		cat("Subst.Freqs = 0.2,0.3,0.3,0.2")
		cat("\n")
		cat("Indel.Model = geo")
		cat("\n")
		cat("Indel.Params =", IM_experimental)
		cat("\n")
		cat("Indel.Rate =", update_ir)
		cat("\n")
		cat("Sim.Reps = 100")
		cat("\n")
		sink()

		EstimateFile<-readLines(file("Final.dawg", "r"))
		system ('dawg - -o R.fasta --label=on', intern = FALSE, input = EstimateFile)
		#system ('perl -i -pe \'s/%r/++$i/eg\' R.fasta', intern = FALSE)
		system ('./NgilaWrapper R.fasta',intern = FALSE)	

		sim_values <- read.table ("sim_values.txt", 
		col.names = c("BL_final", "IM_final", "L_final", "nuisance"))
		attach(sim_values)
		
		check_m <- ((IM_exact - IM_final)^2 / IM_exact)

		system('rm Output.fasta', intern = FALSE)
		system('rm sim_values.txt', intern = FALSE)
		detach(sim_values)
		return (check_m)	
	}

	answer_m <- optimize(M_set, interval_m, tol = 0.001)
	update_m <- answer_m$minimum
	check_m <- answer_m$objective

########################################################################################

	L_set <- function(x) {


		data <- file("Target_Values.txt", "r")

		L_exact <- readLines(data, n=1)
		L_exact <- readLines(data, n=1)
		L_exact <- readLines(data, n=1)
		L_exact <- as.numeric(L_exact)

		L_experimental <- x[1]
		
		#IR_experimental <- L_experimental / (2*SET_BL)
		IR_experimental <- L_experimental / (update_t * 2)	

		sink("Final.dawg")
		cat("Tree.Tree = (C_%r:",update_t,",A_%r:",update_t,");")
		cat("\n")
		cat("Root.Length = 1000")
		cat("\n")
		cat("Subst.Model = jc")
		cat("\n")
		cat("Subst.Freqs = 0.2,0.3,0.3,0.2")
		cat("\n")
		cat("Indel.Model = geo")
		cat("\n")
		cat("Indel.Params =", update_m)
		cat("\n")
		cat("Indel.Rate =", IR_experimental)
		cat("\n")
		cat("Sim.Reps = 100")
		cat("\n")
		sink()

		EstimateFile<-readLines(file("Final.dawg", "r"))
		system ('dawg - -o R.fasta --label=on', intern = FALSE, input = EstimateFile)
		#system ('perl -i -pe \'s/%r/++$i/eg\' R.fasta', intern = FALSE)
		system ('./NgilaWrapper R.fasta',intern = FALSE)	

		sim_values <- read.table ("sim_values.txt", 
		col.names = c("BL_final", "IM_final", "L_final", "nuisance"))
		attach(sim_values)

		check_l <- ((L_exact - L_final)^2 / L_exact)

		system('rm Output.fasta', intern = FALSE)
		system('rm sim_values.txt', intern = FALSE)
		detach(sim_values)
		return (check_l)	
	}

	answer_l <- optimize(L_set, interval_l, tol = 0.001)
	update_l <- answer_l$minimum
	check_l <- answer_l$objective

	cat("======================================================== \n")
	cat("\n")
	cat("cycle number: ", count_cycles, " of ", opt$cycles, "\n")
	cat("\n")
	cat("[t estimate, cost]: [", update_t, ", ", check_t, "] \n")
	cat("[m estimate, cost]: [", update_m, ", ", check_m, "] \n")
	cat("[l estimate, cost]: [", update_l, ", ", check_l, "] \n")
	cat("\n")
	cat("======================================================== \n")

	if (check_t <= opt$Tthreshold & check_m <= opt$Mthreshold & check_l <= opt$Lthreshold) {
		finished <- 1
	}

	count_cycles <- count_cycles + 1
	if (count_cycles > opt$cycles) {
	
		interval_t[1] <- interval_t[1] / 1.5
		interval_t[2] <- interval_t[2] * 1.5

		interval_m[1] <- interval_m[1] / 1.5
		interval_m[2] <- interval_m[2] * 1.5

		interval_l[1] <- interval_l[1] / 1.5
		interval_l[2] <- interval_l[2] * 1.5
		
		opt$cycles <- opt$cycles * 1.5 
		opt$cycles <- round(opt$cycles, 0)
		
		cat("\n")
		cat("the prediction intervals have been widened \n")
		cat("number of cycles has been updated to: ", opt$cycles, "\n")
		cat("\n")
		count_cycles <- 1
	}
}

if (is.null(opt$output) == F) {
	write(c(update_t, update_m, update_ir), file = opt$output,
			ncolumns = 3, append = TRUE, sep = "\t")
}
			
########################################################################################

system('rm Target_Values.txt', intern = FALSE)
system('rm Final.dawg', intern = FALSE)

cat("\n")
cat("======================================================== \n")
cat("\n")
cat("The Simulation is Complete \n")
cat("\n")
cat("Estimated value for Parameter t: ", update_t,"\n") 
cat("Estimated value for Parameter m: ", update_m,"\n") 
cat("Estimated value for Parameter l: ", update_l,"\n") 
cat("\n")
cat("Duration of simulation procedure (s) : \n")
cat("\n")
proc.time() - start_time
cat("\n")
cat("======================================================== \n")
########################################################################################
