function res = contSH_main()

[kp,ic, fp, sd, cp,sp,opt]  = init_par();
[res,nres]                  = init_res(fp,cp,ic,opt);
kp                          = init_idrs(nres,kp,ic,fp);

sd                          = open_log(sd);

res                         = run_cont(res,nres,cp,ic,kp,sd,fp,sp,opt);

                              save_and_plot(res,kp,cp,sp,ic,fp,sd,opt)

end

