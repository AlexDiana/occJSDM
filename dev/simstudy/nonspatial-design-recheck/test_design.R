args<-commandArgs(trailingOnly=TRUE);here<-normalizePath(if(length(args))args[1] else '.')
source(file.path(here,'design_helpers.R'))
root<-normalizePath(file.path(here,'../nonspatial-bias-recheck-20260914'))
load_frozen(root)
a<-readRDS(file.path(root,'operational-k6/qnear_K6-data-01.rds'))
for(arm in c('field4','sites300','knownU')) {
 b<-make_design_input(a,arm)
 check_nested_input(a,b,arm)
 stopifnot(identical(b,make_design_input(a,arm)))
}
far<-readRDS(file.path(root,'operational-k6/qfar_K6-data-01.rds'))
for(arm in c('field4','sites300')) {
 b<-make_design_input(a,arm);c<-make_design_input(far,arm)
 stopifnot(identical(b$sim$true_params$z_true,c$sim$true_params$z_true),identical(b$sim$true_params$w_true,c$sim$true_params$w_true),identical(b$sim$true_params$jsdmParams_true,c$sim$true_params$jsdmParams_true),identical(b$sim$data_list$info,c$sim$data_list$info))
}
# Independent environmental predictions on the original raw scale.
b<-make_design_input(a,'sites300');cols<-paste0('X_psi.EnvCov.',1:2)
ia<-a$sim$data_list$info;ib<-b$sim$data_list$info
xa<-as.matrix(ia[!duplicated(ia$Site),cols]);xb<-as.matrix(ib[!duplicated(ib$Site),cols]);ja<-a$sim$true_params$jsdmParams_true;jb<-b$sim$true_params$jsdmParams_true
raw_beta<-sweep(ja$B,1,apply(xa,2,sd),'/');raw_intercept<-ja$B0-as.vector(colMeans(xa)%*%raw_beta)
independent_eta<-sweep(xb%*%raw_beta+jb$U%*%ja$L,2,raw_intercept,'+')
stopifnot(max(abs(independent_eta-jb$eta))<1e-12)
# Corruption checks must fail before a model can be fitted.
bad<-make_design_input(a,'field4');bad$sim$data_list$OTU[1,1]<-bad$sim$data_list$OTU[1,1]+1
stopifnot(inherits(try(check_nested_input(a,bad,'field4'),silent=TRUE),'try-error'))
ns<-asNamespace('occJSDM');before<-serialize(list(get('runOccJSDM',ns),get('update_jSDMcoef',ns)),NULL)
f<-known_u_fitter(ja$U)
stopifnot(identical(before,serialize(list(get('runOccJSDM',ns),get('update_jSDMcoef',ns)),NULL)))
cat('PASS: nested inputs, deterministic streams, cross-contamination pairing, independent coefficient transformation, corrupt-input rejection, private control leaves installed functions unchanged.\n')
