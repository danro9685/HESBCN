H-ESBCNs (Hidden Extended Suppes-Bayes Causal Networks)
===============================

We implement a C program intended for inferring relationships between cancer mutation events, given an input dataset. Those relationships are depicted in an Extended Suppes-Bayes Causal Network. Furthermore, the method also estimates an HMM to describe survival times. 

With current Makefile, program is compiled with gcc-5. Outputs move information to the console as program runs, with the current sample number output every 1000 iterations. Best poset, along with theta types, printed upon completion of MCMC. Program options are as follows: 

-d | --dataset (required): Path to input file with sample data. Include file type extension as well. ex: input.txt 

-o | --output : Path to desired output file. Final logging information will be written to specified file, as provided by user. If no output specified, results will not be written to file. 

-n | --number_samples : Number of MCMC iterations to run. Default value set to 100,000. Aggressive move threshold and absolute iterations number are factors of 1,000 less/greater than this number. 

-s | --seed : Set seed for random number generator. Default set to time. 

-r | --reg : Set likelihood regularization scheme between "bic", "aic", and "loglik". Default set to "bic". 

-h | --help : Prints help message, quits. 

As indicated, input and output files must be specified. Sample usage: 

"./h-esbcn -d input.txt -o output.txt" -- takes input.txt as input dataset, writes output to output.txt. 

"./h-esbcn -d input.txt -o results/output.txt -n 25000 -s 1" -- takes input.txt as input dataset, writes to output.txt in results directory (inside current directory), sets number of MCMC moves to 25,000, and sets seed to 1. 

"./h-esbcn -d input.txt -o output.txt -r aic" -- takes input.txt as input dataset, writes to output.txt, changes regularization scheme to AIC. 

**NETWORK PARAMETER LEARNING** 

Besides learning the structure of the Suppes-Bayes Causal Network, **H-ESBCN** also runs a nested EM algorithm to estimate exponential distribution parameters for each mutation event and dataset error rate. This takes the form of an alternating maximization algorithm, maximizing lambdas and epsilon (error) given the other. Code for EM stored in markov.h file. 

After learning network structure via metropolis-hastings, we add one event to account for diagnostc event, preprocess the dataset for time computation, and enumerate the valid genotypes given the inferred poset. EM uses a pseudo-warm-start, where error rate assumed to be 0 and lambda values set to MLE given error rate = 0. Our code is based on h-cbn by Gerstung et al. (2007), which was designed for evaluation on conjuncive Bayesian Networks. We also use censored expectation to account for diagnstic censoring, in addition to noise. 

To account for model changes, we move precompute "parent lift" arrays and store the compressed integer representations in the global LIFT array. We use these arrays to determine whether parent sets are fulfilled, should therefore contribute to an event's expected time. 

We resolve deadlock by using the parent lifts to modify the expectation computation and a designated "locked_out" function to test for lockout. Multiple end states are accounted for through genetic completion, and by considering only genotypes with events present into each event's expectation. 

We also use the "is_after" function to determine whether two genotypes lie on the same lineage path, and if so, which comes first. This is integral to proper lineage comparisons. 
