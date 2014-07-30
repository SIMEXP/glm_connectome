function [in,out,opt] = glmc_brick_simple_simus_fdr(in,out,opt)

%% Output
if (nargin < 2)||isempty(out)
    out = 'gb_niak_omitted';
end

%% Options
list_opt = { 'psi' , 'theta' , 'flag_verbose' , 'flag_test' , 'rand_seed' , 'n' , 'm' , 'tx' , 'ty' , 'p_sig'  , 'nb_samp' , 'q'  , 'type_fdr' , 'type_model' };
list_def = { []    , []      , true           , false       , []          , 100 , []  , 0    , 0    , 0.000005 , 1000      , 0.05 , 'GBH'      , 'gaussian'   };
if nargin < 3
    opt = psom_struct_defaults(struct(),list_opt,list_def);
else
    opt = psom_struct_defaults(opt,list_opt,list_def);
end

if isempty(opt.m)
    opt.m = opt.n;
end

if strcmp(opt.type_model,'gaussian')&&isempty(opt.theta) 
    error('Please specify OPT.THETA to use the model ''gaussian''')
end

if strcmp(opt.type_model,'gaussian')&&isempty(opt.theta) 
    error('Please specify OPT.PSI to use the model ''gaussian''')
end

%% If the test flag is true, stop here !
if opt.flag_test == 1
    return
end

%% Seed the random generator
if ~isempty(opt.rand_seed)
    psom_set_rand_seed(opt.rand_seed);
end
nb_q = length(opt.q);
fdr = zeros([opt.nb_samp nb_q]);
sens = zeros([opt.nb_samp nb_q]);
spec = zeros([opt.nb_samp nb_q]);
nb_disc = zeros([opt.nb_samp nb_q]);
for num_q = 1:length(opt.q)
    if opt.flag_verbose
        fprintf('q = %1.3f\n',opt.q(num_q));
    end
    for num_s = 1:opt.nb_samp
        if opt.flag_verbose
            niak_progress(num_s,opt.nb_samp);
        end   
        switch opt.type_model
            case 'pce'
                pce = rand(opt.m,opt.n);
                mask = false(size(pce));
                mask(1:opt.tx,1:opt.ty) = true;
                pce(mask) = opt.p_sig;   
            
            case 'gaussian'
                y = zeros(opt.m,opt.n);
                z0 = randn(1);
                mask = false(size(y));
                for num_g = 1:length(opt.tx)
                    mask1 = false(size(y));
                    mask2 = false(size(y));                
                    if num_g == 1
                        mask1(1:opt.tx(num_g),1:opt.ty(num_g)) = true;            
                        mask2(:,1:opt.ty(num_g)) = true;            
                    else 
                        mask1(1:opt.tx(num_g),(opt.ty(num_g-1)+1):(opt.ty(num_g-1)+opt.ty(num_g))) = true;
                        mask2(:,(opt.ty(num_g-1)+1):(opt.ty(num_g-1)+opt.ty(num_g))) = true;
                    end               
                    mask = mask | mask1;
                    y(mask1) = opt.theta(num_g);
                    y(mask2) = y(mask2) + sqrt(1-opt.psi(num_g)) * randn(opt.m*opt.ty(num_g),1) + sqrt(opt.psi(num_g)) * z0;                
                 end
                 pce = 2*(1-normcdf(abs(y),0,1));
             case 'correlation'
                 y = randn(100,opt.m);
                 z0 = randn(100,1);
                 y(:,1:opt.ty(1)) = sqrt(1-opt.psi(1))*y(:,1:opt.ty(1)) + sqrt(opt.psi(1))*repmat(z0,[1 opt.ty(1)]);
                 y = niak_build_correlation(y,true);
                 y = niak_fisher(y)*sqrt(97);             
                 mask = false(opt.m,opt.m);
                 mask(1:opt.ty(1),1:opt.ty(1)) = true;
                 pce = 2*(1-normcdf(abs(y),0,1));
                 pce = niak_vec2mat(pce);
             otherwise 
                 error('%s is an unknown type of simulation',opt.type_model)
        end
        switch opt.type_fdr
            case 'BH'
                [tmp,tmp2] = niak_fdr(pce(:),opt.type_fdr,opt.q(num_q));
            case {'TST','LSL','TST_sym','LSL_sym'}
                [tmp,tmp2] = niak_glm_fdr(niak_mat2vec(pce),opt.type_fdr,opt.q(num_q));
            otherwise
                [tmp,tmp2] = niak_fdr(pce,opt.type_fdr,opt.q);
        end
        if sum(tmp2(:))>0
            fdr(num_s,num_q)     = sum(tmp2(~mask))/sum(tmp2(:));
            sens(num_s,num_q)    = sum(tmp2(mask))/max(sum(mask(:)),1);
            spec(num_s,num_q)    = 1-sum(tmp2(~mask))/max(sum(~mask(:)),1);
            nb_disc(num_s,num_q) = sum(tmp2(:));
        else 
            fdr(num_s,num_q)     = 0;
            sens(num_s,num_q)    = 0;
            spec(num_s,num_q)    = 1;
            nb_disc(num_s,num_q) = 0;
        end
    end
end

if ~strcmp(out,'gb_niak_omitted')
    q = opt.q;
    save(out,'fdr','sens','spec','nb_disc','q');
end

if opt.flag_verbose
    fprintf('FDR: %1.3f, Sensitivity: %1.3f, Specificity: %1.3f\n',mean(fdr),mean(sens),mean(spec));
end
