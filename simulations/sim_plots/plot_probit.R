rm(list = ls())
library(latex2exp); library(ggpubr)
load("simulations/sim_mcmc_output/out_probit_n50.Rdata")
load("simulations/sim_mcmc_output/out_probit_n100.Rdata")
load("simulations/sim_mcmc_output/out_probit_n300.Rdata")
load("simulations/sim_mcmc_output/out_probit_n500.Rdata")


###########################################################################

source("simulations/sim_plots/plots_template.R")

beta_names = list('1' = TeX("$\\beta_{1}$"), 
                  '2' = TeX("$\\beta_{2} $"), 
                  '3' = TeX("$\\beta_{3} $"))
g_names = list('g' = TeX("$g$"))

beta_labeller <- function(variable,value){
  return(beta_names[value])
}
g_labeller <- function(variable,value){
  return(g_names[value])
}


###########################################################################

p = 3; n = c(50, 100, 300, 500)
r = 100
n_num = length(n)
method = factor(1:4)
# method = c("GP-DP", "BP-DP", "GP-N", "BP-N")
n_meth = length(method)


str(out1)
ind.cov = 16 + 1:4
cov.beta = rbind(t(sapply(ind.cov, function(k) colMeans(out1[[k]][1:r,]))),
                 t(sapply(ind.cov, function(k) colMeans(out2[[k]][1:r,]))),
                 t(sapply(ind.cov, function(k) colMeans(out3[[k]][1:r,]))),
                 t(sapply(ind.cov, function(k) colMeans(out4[[k]][1:r,]))))
# cov.beta = data.frame(cov.beta)

df.cov.prob = data.frame(CP = c(cov.beta), beta = rep(1:3, each = n_meth*n_num),
                         n = rep(n, each = n_meth, times = p), 
                         method = rep(method, times = p*n_num))


g.CP = ggplot(df.cov.prob, aes(x = as.factor(n), y = CP, color = method, fill = method)) + 
  facet_wrap(beta ~ ., labeller=beta_labeller)+
  geom_bar(position="dodge", stat="identity", width = 0.6) +
  geom_hline(aes(yintercept = 0.95), linetype = "dotted") +
  scale_color_manual(values = c( "dodgerblue4", "grey30","orange3", "sienna4"), 
                     labels = c('GP-DP', 'BP-DP', 'GP-N', 'BP-N')) +
  scale_fill_manual(values = c( 'cornflowerblue', "grey","bisque", "#d88546"), 
                    labels = c('GP-DP', 'BP-DP', 'GP-N', 'BP-N')) +
  xlab("n")  + ylab("Coverage\nprobability") + themegg 
g.CP

mse.mat = rbind(out1$MSE.beta[1:r,], out1$MSE.beta.BP[1:r,], 
                out1$MSE.beta.norm[1:r,], out1$MSE.beta.BP.norm[1:r,],
                out2$MSE.beta[1:r,], out2$MSE.beta.BP[1:r,], 
                out2$MSE.beta.norm[1:r,], out2$MSE.beta.BP.norm[1:r,],
                out3$MSE.beta[1:r,], out3$MSE.beta.BP[1:r,], 
                out3$MSE.beta.norm[1:r,], out3$MSE.beta.BP.norm[1:r,],
                out4$MSE.beta[1:r,], out4$MSE.beta.BP[1:r,], 
                out4$MSE.beta.norm[1:r,], out4$MSE.beta.BP.norm[1:r,])


mse.mat = rbind(out1$MSE.beta[1:r,], out1$MSE.beta.BP[1:r,], 
                out1$MSE.beta.norm[1:r,], out1$MSE.beta.BP.norm[1:r,],
                out2$MSE.beta[1:r,], out2$MSE.beta.BP[1:r,], 
                out2$MSE.beta.norm[1:r,], out2$MSE.beta.BP.norm[1:r,],
                out3$MSE.beta[1:r,], out3$MSE.beta.BP[1:r,], 
                out3$MSE.beta.norm[1:r,], out3$MSE.beta.BP.norm[1:r,],
                out4$MSE.beta[1:r,], out4$MSE.beta.BP[1:r,], 
                out4$MSE.beta.norm[1:r,], out4$MSE.beta.BP.norm[1:r,])

df.mse.beta = data.frame(mse = c(mse.mat), 
                         method = rep(method, each = r, times = n_num*p),
                         n = rep(n, each = r*n_meth, times = p), 
                         beta = rep(1:3, each = r*n_meth*n_num))

g.mse.beta = ggplot(df.mse.beta, aes(x = as.factor(n), y = mse, col = method)) + 
  geom_boxplot(outlier.size = 0.2) +
  facet_wrap(beta ~ ., labeller=beta_labeller)+
  scale_color_manual(values = c( "dodgerblue", "grey40","orange", "sienna"), 
                     labels = c('GP-DP', 'BP-DP', 'GP-N', 'BP-N'))+
  xlab("n") + ylab("Mean\nsquared error") + coord_cartesian(ylim = c(0,0.01)) + themegg 
g.mse.beta

dat.mise.g = data.frame(mise = c(out1$MISE.g[1:r]  , out1$MISE.g.BP[1:r] ,
                                 out1$MISE.g.norm[1:r], out1$MISE.g.BP.norm[1:r],
                                 out2$MISE.g[1:r]  , out2$MISE.g.BP[1:r]  ,
                                 out2$MISE.g.norm[1:r], out2$MISE.g.BP.norm[1:r],
                                 out3$MISE.g[1:r]  , out3$MISE.g.BP[1:r]  ,
                                 out3$MISE.g.norm[1:r], out3$MISE.g.BP.norm[1:r],
                                 out4$MISE.g[1:r]  , out4$MISE.g.BP[1:r],
                                 out4$MISE.g.norm[1:r], out4$MISE.g.BP.norm[1:r]) ,
                        method = rep(method, each = r, times = n_num), g = "g",
                        n = rep(n, each = r*n_meth))
g.mise.g = ggplot(dat.mise.g, aes(x = as.factor(n), y = mise, col = method)) + 
  geom_boxplot(outlier.size = 0.2) + facet_wrap(g~., labeller=g_labeller)+
  scale_color_manual(values = c( "dodgerblue", "grey40","orange", "sienna"), 
                     labels = c('GP-DP', 'BP-DP', 'GP-N', 'BP-N'))+ 
  coord_cartesian(ylim = c(0,0.26))+
  xlab("n") + ylab("Mean integrated\nsquared error") + themegg 
g.mise.g


rb.mat = rbind(out1$RB.beta[1:r,], out1$RB.beta.BP[1:r,], 
               out1$RB.beta.norm[1:r,], out1$RB.beta.BP.norm[1:r,],
               out2$RB.beta[1:r,], out2$RB.beta.BP[1:r,], 
               out2$RB.beta.norm[1:r,], out2$RB.beta.BP.norm[1:r,],
               out3$RB.beta[1:r,], out3$RB.beta.BP[1:r,], 
               out3$RB.beta.norm[1:r,], out3$RB.beta.BP.norm[1:r,],
               out4$RB.beta[1:r,], out4$RB.beta.BP[1:r,], 
               out4$RB.beta.norm[1:r,], out4$RB.beta.BP.norm[1:r,])

df.rb.beta = data.frame(rb = c(rb.mat), 
                        method = rep(method, each = r, times = n_num*p),
                        n = rep(n, each = r*n_meth, times = p), 
                        beta = rep(1:3, each = r*n_meth*n_num))

g.rb.beta = ggplot(df.rb.beta, aes(x = as.factor(n), y = rb, col = method)) + 
  geom_boxplot(outlier.size = 0.2) +
  facet_wrap(beta ~ ., labeller=beta_labeller)+
  geom_hline(aes(yintercept = 0), linetype = "dotted") +
  scale_color_manual(values = c( "dodgerblue", "grey40","orange", "sienna"), 
                     labels = c('GP-DP', 'BP-DP', 'GP-N', 'BP-N'))+
  coord_cartesian(ylim = c(-0.5,0.5))+
  xlab("n") + ylab("Relative\nbias") +  themegg 
g.rb.beta

g.all = ggarrange(g.mse.beta, g.CP, g.rb.beta, g.mise.g, nrow = 4, ncol = 1, common.legend = TRUE)
g.all
ggsave("simulations/sim_plots/figures/simprobit.pdf", g.all, width = 7, height = 9)


