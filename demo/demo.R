## Throughout this, you can see the flag
## "TRY ME!" highlighting areas to play with-- use
## command-F (or windows equivalent) to find!

####################################################
########## Demo Setup --- Run this once! ###########
####################################################
                                                  ##
install.packages(c("remotes", "Rcpp", "abind",    ##
                   "pracma", "sparseEigen",       ##
                   "mvtnorm", "gridExtra",        ##
                   "ggplot2"))                    ##
                                                  ##
remotes::install_github("kelrenmor/bs3fa",        ##
                        quiet=T, upgrade=T,       ##
                        dependencies=T)           ##
                                                  ##
####################################################
####################################################

########################################################
### Demonstration of BS3FA in Dose-Response Modeling ###
########################################################

# --- SLIDE 2 --- #

## Load required libraries -- installed in demoSetup.R
library(ggplot2)
library(gridExtra)
library(pracma)

## Load bs3fa package and dataset drdat
library(bs3fa)
data('drdat')

## Set up values for data simulation

# --- SLIDE 3 --- #

##############
### TRY ME ### Set a different seed (any integer value) for different results
##############
set.seed(12345) # set seed for random number generation

N=300           # no. of 'chemicals'
D=20; S=40      # no. of unique doses (D) and features (S)
K=3; J=4        # dim of shared (K) and X-specific (J) latent space
std_error_y=0.2 # noise in observing curves (Y)
std_error_x=0.1 # noise in observing features (X)

############### Make prob_miss larger/smaller in [0,1):
### TRY ME! ### - See how much missingness it takes to break the model.
############### - See how the uncertainty changes as data are more/less observed.
prob_miss=0.7   # expected proportion of dose-response grid unobserved

real_Y=drdat    # use ToxCast data to simulate basis for Y 

###############
### TRY ME! ### Uncomment the line below to fully simulate data rather
############### than using true chemical data to simulate a basis Lambda.
# real_Y=NULL

# --- SLIDE 4 --- #

## Simulate dose-response and feature data
dat=bs3fa::simulate_data(real_Y=real_Y,N=N,D=D,S=S,K=K,J=J,
                         prob_miss=prob_miss,
                         std_error_y=std_error_y,
                         std_error_x=std_error_x)

## Plot average dose-response curve
qplot(dat$doses, dat$avg_dose_resp, geom=c("point", "line")) + 
  theme_minimal() + theme_bw() + xlab("dose") + ylab("average response")

# --- SLIDE 5 --- #

## Plot example dose-response curves
bs3fa::plot_data(Y=dat$Y, true_curve=dat$Lambda_true%*%dat$eta_true, 
                 avg_dose_resp=dat$avg_dose_resp, doses=dat$doses, inds=1:3)

# --- SLIDE 6 --- #

## Plot true $Y$ factor loadings $\Lambda$
bs3fa::plot_matrix(dat$Lambda_true, type="Lambda")                # (visual 1)
bs3fa::plot_Lambda_true(Lambda=dat$Lambda_true, doses=dat$doses)  # (visual 2)

# --- SLIDE 8 --- #

## Plot true $X$ common factor loadings $\Theta$ 
bs3fa::plot_matrix(dat$Theta_true, type="Theta")                         # (visual 1)
qplot(c(dat$Theta_true), bins=20) + xlab(expression(theta[s*","~i])) +   # (visual 2)
  ylab("frequency") + theme_minimal() + theme_bw()

# --- SLIDE 10 --- #

## Prior to running the model, randomly set 20% of dose-response curves to NA
prop_unobserved = 0.2
not_obs = sample(1:N, round(prop_unobserved*N))
dat$Y[,not_obs] = NA

# --- SLIDE 11 --- #

## Set sampler settings
K_p=dat$K+3 # give sampler 'guess' at true shared loadings
J_p=dat$J+3 # give sampler 'guess' at true X-specific loadings
thin=10     # keep every 10th sample
burnin=5000 # run sampler for 5000 draws before beginning to save
nsamps_save=500 # save 500 samples in total
post_process=T # resolve label/sign switches, rotational ambiguity

## Run Gibbs sampler
res=bs3fa::run_bs3fa(X=dat$X, Y=dat$Y, K=K_p, J=J_p, 
                     thin=thin, nsamps_save=nsamps_save, 
                     burnin=burnin, post_process=post_process,
                     print_progress=T)

# --- SLIDE 16 --- #

## Plot model performance ($\Lambda$) pre-cleaning
## Note how the columns aren't necessarily matched up to truth, and the signs don't match
p1 = bs3fa::plot_matrix(dat$Lambda_true, type="Lambda", tit="Truth", include_legend=F)
Lam_tmp = matrix(NA,nrow=nrow(dat$Lambda_true),ncol=K_p)
for(k_mod in 1:K_p){Lam_tmp[,k_mod] = apply(res$Lambda_save[,k_mod,],1,mean)}
p2 <- bs3fa::plot_matrix(Lam_tmp, type="Lambda", tit="Estimate", include_legend=F)
gridExtra::grid.arrange(p1, p2, nrow = 1, widths = c(3/9,6/9))

# --- SLIDE 17 --- #

## Match model column indices and sign to truth
# When 'truth' is known, re-index model runs to match true indices
# (this is just for later plotting convenience)
res_clean = bs3fa::reorder_entries(eta_true=dat$eta_true, eta=res$eta_save, 
                                   Lambda_true=dat$Lambda_true, Lambda=res$Lambda_save,
                                   Theta_true=dat$Theta_true,Theta=res$Theta_save)

# --- SLIDE 18 --- #

## Plot model performance ($\Theta$) post-cleaning
# Compare truth and model mean side-by-side
p1 <- bs3fa::plot_matrix(dat$Theta_true, type="Theta", tit="Truth", include_legend=F)
p2 <- bs3fa::plot_matrix(res_clean$Theta, type="Theta", tit="Estimate", include_legend=F)
gridExtra::grid.arrange(p1, p2, nrow = 1, widths = c(3/9,6/9))
# Plot lower and upper end of 95% credible interval side-by-side
p1 <- bs3fa::plot_matrix(res_clean$Theta_low, type="Theta", tit="Lower 2.5%", include_legend=F)
p2 <- bs3fa::plot_matrix(res_clean$Theta_upp, type="Theta", tit="Upper 97.5%", include_legend=F)
gridExtra::grid.arrange(p1, p2, nrow = 1)

# --- SLIDE 20 --- #

## Plot model performance ($\Lambda$) post-cleaning
# Compare truth and model mean side-by-side
p1 <- bs3fa::plot_matrix(dat$Lambda_true, type="Lambda", tit="Truth", include_legend=F)
p2 <- bs3fa::plot_matrix(res_clean$Lambda, type="Lambda", tit="Estimate", include_legend=F)
gridExtra::grid.arrange(p1, p2, nrow = 1, widths = c(3/9,6/9))
# Plot lower and upper end of 95% credible interval side-by-side
p1 <- bs3fa::plot_matrix(res_clean$Lambda_low, type="Lambda", tit="Lower 2.5%", include_legend=F)
p2 <- bs3fa::plot_matrix(res_clean$Lambda_upp, type="Lambda", tit="Upper 97.5%", include_legend=F)
gridExtra::grid.arrange(p1, p2, nrow = 1)
# Plot truth and 95% credible interval by dose for indices 1:3
bs3fa::plot_Lambda_mod_tru(Lambda_low=res_clean$Lambda_low, Lambda=res_clean$Lambda,
Lambda_upp=res_clean$Lambda_upp,
Lambda_true=dat$Lambda_true,
doses=dat$doses, inds=1:3)

# --- SLIDE 24 --- #

## Plot model performance ($\eta$)
qplot(c(dat$eta_true), c(res_clean$eta[1:K,]), geom="point", alpha=0.2, show.legend = F) +
xlab(expression(True~eta)) + ylab(expression(Predicted~eta)) + 
geom_abline(intercept=0, slope=1, color="red") + 
theme_minimal() + theme_bw()

# --- SLIDE 25 --- #

## Plot model performance (distance)
dist_tru = c(dist(t(dat$eta_true)))
dist_mod = c(dist(t(res_clean$eta[1:K,])))
samp_inds = sample(1:length(dist_tru), 5*1e2) # sample for plotting purposes
qplot(dist_tru, dist_mod, geom="point", alpha=0.8, show.legend = F) +
xlab("True distance") + ylab(expression("Predicted distance")) + 
geom_abline(intercept=0, slope=1, color="red") + 
theme_minimal() + theme_bw()

# --- SLIDE 26 --- #

## Plot model performance (prediction)
# Turn Lambda and eta samples into Lambda * eta 
# (i.e., the predicted dose response curve)
res_pred = pred_drcurve(Lambda_mod=res$Lambda_save, eta_mod=res$eta_save, rescale=1)
# Display three random example unobserved curves
ind_samp = sample(1:length(not_obs),2) # Pick two random unobserved curves

###############
### TRY ME! ### Uncomment the line below to pick some other curve predictions to view,
############### setting it to any pair of numbers from 1,...,60 (length(not_obs)).
# ind_samp = c(2,40)

miss_ex = not_obs[ind_samp] 
bs3fa::plot_Lambda_mod_tru(Lambda_low=res_pred$dr_low, Lambda=res_pred$dr_est,
                           Lambda_upp=res_pred$dr_upp,Lambda_true=dat$true_curve,
                           doses=dat$doses, inds=miss_ex, ylab_head="response", 
                           title_head="Y", mean_curve=dat$avg_dose_resp)

# --- BONUS CONTENT --- #

## View contributions of individual factors and 
## how they become the true dose-response curve.
i = 1 # Choose chemical index to look at (i=1 for slide 7 example).
dat$eta_true[,i] # The weights associated with each column of Lambda.
piece_1 = dat$Lambda_true[,1] * dat$eta_true[1,i]
piece_2 = dat$Lambda_true[,2] * dat$eta_true[2,i]
piece_3 = dat$Lambda_true[,3] * dat$eta_true[3,i]
plot(dat$doses, piece_1, type="l", xlab="dose", ylab="", main="Response contribution by 1st factor")
plot(dat$doses, piece_2, type="l", xlab="dose", ylab="", main="Response contribution by 2nd factor")
plot(dat$doses, piece_3, type="l", xlab="dose", ylab="", main="Response contribution by 3rd factor")
plot(dat$doses, piece_1 + piece_2 + piece_3 + dat$avg_dose_resp, 
     type="l", xlab="dose", ylab="", main="True curve")

## Get the true and estimated std deviation of Y and X
# For Y
mod_sd_y = sqrt(res$sigsq_y_save)
hist(mod_sd_y)
abline(v=mean(mod_sd_y), col="blue") # model estimate
abline(v=std_error_y, col="red") # truth 
# For X
mod_sd_x = sqrt(res$sigsq_x_save)
hist(mod_sd_x)
abline(v=mean(mod_sd_x), col="blue") # model estimate
abline(v=std_error_x, col="red") # truth 
