#' sr_plot function
#'
#' This function generates a stock-recruitment (S-R) model for rstan based on model inputs.
#' @param type Specify whether to generate a 'static' S-R model, where parameters are time-invariant, 
#' a time-varying 'tv' model, or a regime shift model 'regime'
#' @param df stock-recruitment dataset. With columns R (recruits) and S (spawners) and by (brood year)
#' @param form Either fit with Stan ('stan') or TMB ('tmb')
#' @param mod a fitted Stan or TMB model
#' @param pdf whether to create a pdf or not
#' @return returns the compiled rstan code for a given S-R model
#' @importFrom rstan stan_model
#' @export
#' @examples
#' sr_plot(type='static',df=df,form='stan',df=df,mod=f1,pdf=FALSE)

sr_plot=function(df,mod,title,make.pdf=FALSE,path,type=c('static','rw','hmm'),par=c('a','b','both'),form=c('stan','tmb'),ac=FALSE,sr_only=FALSE){
  if(type=='static'){ #static====
    x_new=seq(0,max(df$S),length.out=200)

    if(form=='stan'){
      post=rstan::extract(mod)
      pred_df=data.frame(pred=exp(median(post$log_a)-median(post$b)*x_new)*x_new,x_new=x_new)
      }
    if(form=='tmb'){
      pred_df=data.frame(pred=exp(mod$alpha-mod$beta*x_new)*x_new,x_new=x_new)
    }
    plot=ggplot2::ggplot(df, aes(S, R)) +
      geom_line(data=pred_df,aes(x=x_new,y=pred),linewidth=1.3)+
      geom_point(aes(colour = by),size=3.5) +
      scale_colour_viridis_c(name='Year')+
      ggtitle(title)+
      xlab("Spawners") + 
      ylab("Recruits")+
      xlim(0, max(df$S))+
      ylim(0, max(df$R))+
      theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
            strip.text = element_text(face="bold", size=12),
            axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
  }
  if(type=='rw'){
    x_new=seq(0,max(df$S),length.out=200)
    by_q=round(quantile(df$by,seq(0,1,by=0.1)))
      if(par=='a'){ #rw alpha=====
        if(form=='stan'){
          post=rstan::extract(mod)
          pred_df=data.frame(x_new)
        for(n in 1:length(by_q)){
          pred_df[,1+n]=exp(median(post$log_a[,match(by_q[n],df$by)])-median(post$b)*x_new)*x_new
        }
          alpha_df=data.frame(by=seq(min(df$by),max(df$by)),med=apply(post$log_a,2,median),l90=apply(post$log_a,2,quantile,0.1),u90=apply(post$log_a,2,quantile,0.9))
          plot2=ggplot2::ggplot(alpha_df, aes(by,med)) +
            geom_line(aes(x=by,y=med),linewidth=1.3)+
            geom_point(aes(colour = by),size=4) +
            scale_colour_viridis_c(name='Year')+
            geom_ribbon(aes(ymin =l90, ymax =u90), alpha = 0.2)+
            xlab("Year") + 
            ylab("log(Alpha)")+
            theme_classic(14)+
            theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                  strip.text = element_text(face="bold", size=12),
                  axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          
        }
        
        if(form=='tmb'){
          pred_df=data.frame(x_new)
          for(n in 1:length(by_q)){
            pred_df[,1+n]=exp(mod$alpha[match(by_q[n],df$by)]-mod$beta*x_new)*x_new
          }
          alpha_df=data.frame(by=df$by,med=mod$alpha)
          
          plot2=ggplot2::ggplot(alpha_df, aes(by,med)) +
            geom_hline(yintercept=0,linetype='dashed')+
            geom_line(aes(x=by,y=med),linewidth=1.3)+
            geom_point(aes(colour = by),size=4) +
            scale_colour_viridis_c(name='Year')+
            xlab("Year") + 
            ylab(paste0("Productivity - log(","\u03b1",")"))+
            theme_classic(14)+
            theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                  strip.text = element_text(face="bold", size=12),
                  axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          
          
        }
        
        plot1=ggplot2::ggplot(df, aes(S/1000, R/1000)) +
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = by_q[1]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = by_q[2]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,4],colour = by_q[3]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,5],colour = by_q[4]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,6],colour = by_q[5]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,7],colour = by_q[6]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = by_q[1]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = by_q[2]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,4],colour = by_q[3]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,5],colour = by_q[4]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,6],colour = by_q[5]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,7],colour = by_q[6]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,8],colour = by_q[7]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,9],colour = by_q[8]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,10],colour = by_q[9]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,11],colour = by_q[10]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,12],colour = by_q[11]),linewidth=1.3)+
          geom_point(aes(colour = by),size=2.5) +
          scale_colour_viridis_c(name='Year')+
          xlab("Spawners (thousands)") + 
          ylab("Recruits (thousands)")+
          xlim(0, max(df$S))+
          ylim(0, max(df$R))+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))

       
        title <- ggdraw() + 
          draw_label(
            title,
            fontface = 'bold',
            x = 0,
            hjust = 0
          ) +
          theme(
            # add margin on the left of the drawing canvas,
            # so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 0, 7)
          )
        
        legend = cowplot::get_legend(plot1)
        
        plot_rw_a=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                      plot2 + theme(legend.position="none"),
                       ncol=2,nrow=1)
        plot=cowplot::plot_grid(plot_rw_a,legend,rel_widths = c(3,.25))
        plot=cowplot::plot_grid(title,plot,ncol=1,rel_heights = c(0.1,1))
        if(sr_only==TRUE){plot=plot1}
      }
      if(par=='b'){ ###rw beta=====
          if(form=='stan'){
            post=rstan::extract(mod)
            pred_df=data.frame(x_new)
            
            for(n in 1:length(by_q)){
              pred_df[,1+n]=exp(median(post$log_a)-median(post$b[,match(by_q[n],df$by)])*x_new)*x_new
            }
            
            beta_df=data.frame(by=seq(min(df$by),max(df$by)),med=apply(post$S_max,2,median),l90=apply(post$S_max,2,quantile,0.1),u90=apply(post$S_max,2,quantile,0.9))
            
            plot2=ggplot2::ggplot(beta_df, aes(by,med)) +
              geom_line(aes(x=by,y=med),linewidth=1.3)+
              geom_point(aes(colour = by),size=4) +
              scale_colour_viridis_c(name='Year')+
              geom_ribbon(aes(ymin =l90, ymax =u90), alpha = 0.2)+
              ggtitle(paste(title))+
              xlab("Year") + 
              ylab("Smax")+
              theme_classic(14)+
              theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                    strip.text = element_text(face="bold", size=12),
                    axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
            
          }
          if(form=='tmb'){
            pred_df=data.frame(x_new)
            for(n in 1:length(by_q)){
              pred_df[,1+n]=exp(mod$alpha-mod$beta[match(by_q[n],df$by)]*x_new)*x_new
            }
            beta_df=data.frame(by=df$by,med=mod$Smax)
            
            plot2=ggplot2::ggplot(beta_df, aes(by,med)) +
             geom_line(aes(x=by,y=med),linewidth=1.3)+
              geom_point(aes(colour = by),size=4) +
              scale_colour_viridis_c(name='Year')+
              ggtitle(paste(title))+
              xlab("Year") + 
              ylab("Smax")+
              theme_classic(14)+
              theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                    strip.text = element_text(face="bold", size=12),
                    axis.text=element_text(face="bold",size=14),axis.title = element_text(face="bold",size=14),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
            
            
          }
          
        
        plot1=ggplot2::ggplot(df, aes(S, R)) +
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = by_q[1]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = by_q[2]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,4],colour = by_q[3]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,5],colour = by_q[4]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,6],colour = by_q[5]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,7],colour = by_q[6]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = by_q[1]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = by_q[2]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,4],colour = by_q[3]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,5],colour = by_q[4]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,6],colour = by_q[5]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,7],colour = by_q[6]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,8],colour = by_q[7]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,9],colour = by_q[8]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,10],colour = by_q[9]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,11],colour = by_q[10]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,12],colour = by_q[11]),linewidth=1.3)+
          geom_point(aes(colour = by),size=2.5) +
          scale_colour_viridis_c(name='Year')+
          ggtitle(title)+
          xlab("Spawners") + 
          ylab("Recruits")+
          xlim(0, max(df$S))+
          ylim(0, max(df$R))+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold",size=14),axis.title = element_text(face="bold",size=14),plot.title = element_text(face = "bold", hjust = 0.5,size=15))

        
        legend = cowplot::get_legend(plot1)
        
        plot_rw_b=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                                     plot2 + theme(legend.position="none"),
                                     ncol=2,nrow=1,labels=c("A","B"))
        plot=cowplot::plot_grid(plot_rw_b,legend,rel_widths = c(3,.25))
        if(sr_only==TRUE){plot=plot1}
      }
      if(par=='both'){ #rw alpha beta=====
        if(form=='stan'){
          post=rstan::extract(mod)
          pred_df=data.frame(x_new)
          for(n in 1:length(by_q)){
            pred_df[,1+n]=exp(median(post$log_a[,match(by_q[n],df$by)])-median(post$b[,match(by_q[n],df$by)])*x_new)*x_new
          }
          
          alphabeta_df=data.frame(by=seq(min(df$by),max(df$by)),a_med=apply(post$log_a,2,median),a_l90=apply(post$log_a,2,quantile,0.15),a_u90=apply(post$log_a,2,quantile,0.85),b_med=apply(post$S_max,2,median),b_l90=apply(post$S_max,2,quantile,0.1),b_u90=apply(post$S_max,2,quantile,0.9))
         
          plot2=ggplot2::ggplot(alphabeta_df, aes(by,a_med)) +
            geom_line(aes(x=by,y=a_med),linewidth=1.3)+
            geom_point(aes(colour = by),size=3) +
            scale_colour_viridis_c(name='Year')+
            geom_ribbon(aes(ymin =a_l90, ymax =a_u90), alpha = 0.2)+
            ggtitle(paste(title))+
            xlab("Year") + 
            ylab("log(Alpha)")+
            theme_classic(14)+
            theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                  strip.text = element_text(face="bold", size=12),
                  axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          
          plot3=ggplot2::ggplot(alphabeta_df, aes(by,b_med)) +
            geom_line(aes(x=by,y=b_med),linewidth=1.3)+
            geom_point(aes(colour = by),size=3) +
            scale_colour_viridis_c(name='Year')+
            geom_ribbon(aes(ymin =b_l90, ymax =b_u90), alpha = 0.2)+
            ggtitle(paste(title))+
            xlab("Year") + 
            ylab("Smax")+
            theme_classic(14)+
            theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                  strip.text = element_text(face="bold", size=12),
                  axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          
        }
        if(form=='tmb'){
          pred_df=data.frame(x_new)
          for(n in 1:length(by_q)){
            pred_df[,1+n]=exp(mod$alpha[match(by_q[n],df$by)]-mod$beta[match(by_q[n],df$by)]*x_new)*x_new
          }
          alphabeta_df=data.frame(by=df$by,a_med=mod$alpha,b_med=mod$Smax)
          
          plot2=ggplot2::ggplot(alphabeta_df, aes(by,a_med)) +
            geom_line(aes(x=by,y=a_med),linewidth=1.3)+
            geom_point(aes(colour = by),size=3) +
            scale_colour_viridis_c(name='Year')+
            ggtitle(paste(title))+
            xlab("Year") + 
            ylab("log(Alpha)")+
            theme_classic(14)+
            theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                  strip.text = element_text(face="bold", size=12),
                  axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          
          plot3=ggplot2::ggplot(alphabeta_df, aes(by,b_med)) +
            geom_line(aes(x=by,y=b_med),linewidth=1.3)+
            geom_point(aes(colour = by),size=3) +
            scale_colour_viridis_c(name='Year')+
            ggtitle(paste(title))+
            xlab("Year") + 
            ylab("Smax")+
            theme_classic(14)+
            theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                  strip.text = element_text(face="bold", size=12),
                  axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          
          
        }
        
        plot1=ggplot2::ggplot(df, aes(S, R)) +
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = by_q[1]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = by_q[2]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,4],colour = by_q[3]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,5],colour = by_q[4]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,6],colour = by_q[5]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,7],colour = by_q[6]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = by_q[1]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = by_q[2]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,4],colour = by_q[3]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,5],colour = by_q[4]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,6],colour = by_q[5]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,7],colour = by_q[6]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,8],colour = by_q[7]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,9],colour = by_q[8]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,10],colour = by_q[9]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,11],colour = by_q[10]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,12],colour = by_q[11]),linewidth=1.3)+
          geom_point(aes(colour = by),size=2.5) +
          scale_colour_viridis_c(name='Year')+
          ggtitle(title)+
          xlab("Spawners") + 
          ylab("Recruits")+
          xlim(0, max(df$S))+
          ylim(0, max(df$R))+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
        
        
        legend = cowplot::get_legend(plot1)
        
        plot=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                                     plot2 + theme(legend.position="none"),
                                     plot3 + theme(legend.position="none"),
                                     legend,
                                     ncol=2,nrow=2,labels=c("A","B","C"))
        if(sr_only==TRUE){plot=plot1}
        
      }
      
    }
    if(type=='hmm'){ 
        
        x_new=seq(min(df$S),max(df$S),length.out=200)
        pred_df=data.frame(x_new)
        
        if(par=='a'){ #hmm alpha====
          
          if(form=='stan'){
          post=rstan::extract(mod)
          pred_df[,2]=exp(median(post$log_a[,1])-median(post$b)*x_new)*x_new
          pred_df[,3]=exp(median(post$log_a[,2])-median(post$b)*x_new)*x_new
          df$gamma=apply(post$gamma[,,2],2,median)
          gamma_df=data.frame(by=df$by,gamma=apply(post$gamma[,,2],2,median),gamma_l90=apply(post$gamma[,,2],2,quantile,0.1),gamma_u90=apply(post$gamma[,,2],2,quantile,0.9))
          
          plot2=ggplot2::ggplot(gamma_df, aes(by,gamma)) +
            ylim(0,1)+
            geom_hline(yintercept=0.5,linetype='dashed')+
            geom_line(aes(x=by,y=gamma),linewidth=1.3)+
            geom_point(aes(colour = gamma),size=4) +
            scale_colour_viridis_c(name='p')+
            geom_ribbon(aes(ymin =gamma_l90, ymax =gamma_u90), alpha = 0.2)+
            ggtitle(paste(title))+
            xlab("Year") + 
            ylab("Prob. of high prod. regime")+
            theme_classic(14)+
            theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                  strip.text = element_text(face="bold", size=12),
                  axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          }
          
          if(form=='tmb'){
            pred_df[,2]=exp(mod$alpha[1]-mod$beta*x_new)*x_new
            pred_df[,3]=exp(mod$alpha[2]-mod$beta*x_new)*x_new
            gamma_df=data.frame(by=df$by,gamma=mod$probregime[2,])
  
            plot2=ggplot2::ggplot(gamma_df, aes(by,gamma)) +
              ylim(0,1)+
              geom_hline(yintercept=0.5,linetype='dashed')+
              geom_line(aes(x=by,y=gamma),linewidth=1.3)+
              geom_point(aes(colour = gamma),size=4) +
              scale_colour_viridis_c(name='P(High prod. regime)')+
              ggtitle(paste(title))+
              xlab("Year") + 
              ylab("Prob. of high prod. regime")+
              theme_classic(14)+
              theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                    strip.text = element_text(face="bold", size=12),
                    axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
            
          }
          
          plot1=ggplot2::ggplot(df, aes(S, R)) +
            geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = min(gamma_df$gamma)),linewidth=1.3)+
            geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = max(gamma_df$gamma)),linewidth=1.3)+
            geom_point(aes(colour = gamma_df$gamma),size=2.5) +
            scale_colour_viridis_c(name='p')+
            ggtitle(title)+
            xlab("Spawners") + 
            ylab("Recruits")+
            theme_classic(14)+
            xlim(0, max(df$S))+
            ylim(0, max(df$R))+
            theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                  strip.text = element_text(face="bold", size=12),
                  axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          
          
          
          legend = cowplot::get_legend(plot1)
          
          plot_hmm_a=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                                       plot2 + theme(legend.position="none"),
                                       ncol=2,nrow=1,labels=c("A","B"))
          plot=cowplot::plot_grid(plot_hmm_a,legend,rel_widths = c(3,.3))
          if(sr_only==TRUE){plot=plot1}
        }
        
        if(par=='b'){ #hmm beta====
          
          if(form=='stan'){
          post=rstan::extract(mod)
          pred_df[,2]=exp(median(post$log_a)-median(post$b[,1])*x_new)*x_new
          pred_df[,3]=exp(median(post$log_a)-median(post$b[,2])*x_new)*x_new
          df$gamma=apply(post$gamma[,,2],2,median)
          gamma_df=data.frame(by=df$by,gamma=apply(post$gamma[,,1],2,median),gamma_l90=apply(post$gamma[,,1],2,quantile,0.1),gamma_u90=apply(post$gamma[,,1],2,quantile,0.9))
          
          plot2=ggplot2::ggplot(gamma_df, aes(by,gamma)) +
            ylim(0,1)+
            geom_hline(yintercept=0.5,linetype='dashed')+
            geom_line(aes(x=by,y=gamma),linewidth=1.3)+
            geom_point(aes(colour = gamma),size=4) +
            scale_colour_viridis_c(name='p')+
            geom_ribbon(aes(ymin =gamma_l90, ymax =gamma_u90), alpha = 0.2)+
            ggtitle(paste(title))+
            xlab("Year") + 
            ylab("Prob. of high cap. regime")+
            theme_classic(14)+
            theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                  strip.text = element_text(face="bold", size=12),
                  axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          }
          if(form=='tmb'){
            pred_df[,2]=exp(mod$alpha-mod$beta[1]*x_new)*x_new
            pred_df[,3]=exp(mod$alpha-mod$beta[2]*x_new)*x_new
            gamma_df=data.frame(by=df$by,gamma=mod$probregime[1,])
            
            plot2=ggplot2::ggplot(gamma_df, aes(by,gamma)) +
              ylim(0,1)+
              geom_hline(yintercept=0.5,linetype='dashed')+
              geom_line(aes(x=by,y=gamma),linewidth=1.3)+
              geom_point(aes(colour = gamma),size=4) +
              scale_colour_viridis_c(name='p')+
              ggtitle(paste(title))+
              xlab("Year") + 
              ylab("Prob. of high cap. regime")+
              theme_classic(14)+
              theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                    strip.text = element_text(face="bold", size=12),
                    axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
            
          }
          
          plot1=ggplot2::ggplot(df, aes(S, R)) +
            geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = max(gamma_df$gamma)),linewidth=1.3)+
            geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = min(gamma_df$gamma)),linewidth=1.3)+
            geom_point(aes(colour = gamma_df$gamma),size=2.5) +
            scale_colour_viridis_c(name='p')+
            ggtitle(title)+
            xlab("Spawners") + 
            ylab("Recruits")+
            xlim(0, max(df$S))+
            ylim(0, max(df$R))+
            theme_classic(14)+
            theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                  strip.text = element_text(face="bold", size=12),
                  axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          
         
          
          legend = cowplot::get_legend(plot1)
          
          plot_hmm_b=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                                        plot2 + theme(legend.position="none"),
                                        ncol=2,nrow=1,labels=c("A","B"))
          plot=cowplot::plot_grid(plot_hmm_b,legend,rel_widths = c(3,.3))
          if(sr_only==TRUE){plot=plot1}
        }
        if(par=='both'){ #hmm alpha beta====
          if(form=='stan'){
          post=rstan::extract(mod)
          pred_df[,2]=exp(median(post$log_a[,1])-median(post$b[,1])*x_new)*x_new
          pred_df[,3]=exp(median(post$log_a[,2])-median(post$b[,2])*x_new)*x_new
          df$gamma=apply(post$gamma[,,2],2,median)
          gamma_df=data.frame(by=df$by,gamma=apply(post$gamma[,,2],2,median),gamma_l90=apply(post$gamma[,,2],2,quantile,0.1),gamma_u90=apply(post$gamma[,,2],2,quantile,0.9))
          
          plot2=ggplot2::ggplot(gamma_df, aes(by,gamma)) +
            ylim(0,1)+
            geom_hline(yintercept=0.5,linetype='dashed')+
            geom_line(aes(x=by,y=gamma),linewidth=1.3)+
            geom_point(aes(colour = gamma),size=4) +
            scale_colour_viridis_c(name='p')+
            geom_ribbon(aes(ymin =gamma_l90, ymax =gamma_u90), alpha = 0.2)+
            ggtitle(paste(title))+
            xlab("Year") + 
            ylab("Prob. of high prod. regime")+
            theme_classic(14)+
            theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                  strip.text = element_text(face="bold", size=12),
                  axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          }
          if(form=='tmb'){
            pred_df[,2]=exp(mod$alpha[1]-mod$beta[1]*x_new)*x_new
            pred_df[,3]=exp(mod$alpha[2]-mod$beta[2]*x_new)*x_new
            gamma_df=data.frame(by=df$by,gamma=mod$probregime[1,])
            
            plot2=ggplot2::ggplot(gamma_df, aes(by,gamma)) +
              ylim(0,1)+
              geom_hline(yintercept=0.5,linetype='dashed')+
              geom_line(aes(x=by,y=gamma),linewidth=1.3)+
              geom_point(aes(colour = gamma),size=4) +
              scale_colour_viridis_c(name='p')+
              ggtitle(paste(title))+
              xlab("Year") + 
              ylab("Prob. of high cap. regime")+
              theme_classic(14)+
              theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                    strip.text = element_text(face="bold", size=12),
                    axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          }
          
          plot1=ggplot2::ggplot(df, aes(S, R)) +
            geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = min(gamma_df$gamma)),linewidth=1.3)+
            geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = max(gamma_df$gamma)),linewidth=1.3)+
            geom_point(aes(colour = gamma_df$gamma),size=2.5) +
            scale_colour_viridis_c(name='p')+
            ggtitle(title)+
            xlab("Spawners") + 
            ylab("Recruits")+
            xlim(0, max(df$S))+
            ylim(0, max(df$R))+
            theme_classic(14)+
            theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                  strip.text = element_text(face="bold", size=12),
                  axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          
          legend = cowplot::get_legend(plot1)
          
          plot_hmm_ab=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                                        plot2 + theme(legend.position="none"),
                                        ncol=2,nrow=1,labels=c("A","B"))
          plot=cowplot::plot_grid(plot_hmm_ab,legend,rel_widths = c(3,.3))
          if(sr_only==TRUE){plot=plot1}
        }
    }
    if(make.pdf==TRUE){
      if(type=='static'&ac==FALSE){ pdf(here(path,paste(paste(title,type,form,sep='_'),'.pdf',sep='')),width=8,height=6)}
      if(type=='static'&ac==TRUE){ pdf(here(path,paste(paste(title,type,'ac',form,sep='_'),'.pdf',sep='')),width=8,height=6)}
      if(type=='rw'&par=='both'){ pdf(here(path,paste(paste(title,type,par,form,sep='_'),'.pdf',sep='')),width=10,height=10)}
      if(type=='rw'&par!='both'){ pdf(here(path,paste(paste(title,type,par,form,sep='_'),'.pdf',sep='')),width=14,height=6)}
      if(type=='hmm'){ pdf(here(path,paste(paste(title,type,par,form,sep='_'),'.pdf',sep='')),width=14,height=6)}
      return(plot)
      
      dev.off()
      dev.off()
    }
  if(make.pdf==FALSE){
    return(plot)
  }
  
  }

#' Prior predictive check 
#'
#' Any arguments to pass to [rstan::sampling].
#'
#' @param fit A model fit in rstan from e.g. ricker_stan, ricker_rw_stan, ricker_hmm_stan
#'
#' @export
prior_check<- function(data,type=c('static','rw'),AC=FALSE,par=c('a','b','both'),control = stancontrol()){
  if(type=='static'&AC==FALSE){
   
    sm=sr_prior_mod(type='static')
      
   if(is.null(smax_priors)==T){
      datm = list(N=nrow(data),
                  L=max(data$by)-min(data$by)+1,
                  ii=seq_len(nrow(data)),
                  R_S =data$logRS,
                  S=data$S,
                  pSmax_mean=max(data$S)/2,
                  pSmax_sig=max(data$S))
      
    }
    if(is.null(smax_priors)==F){
      datm = list(N=nrow(data),
                  L=max(data$by)-min(data$by)+1,
                  ii=seq_len(nrow(data)),
                  R_S =data$logRS,
                  S=data$S,
                  pSmax_mean=smax_priors[1],
                  pSmax_sig=smax_priors[2])
    }
  }
  if(type=='static'&AC==TRUE){
    if(is.null(smax_priors)==T){
      sm=sr_prior_mod(type='static',AC=TRUE)
    }
    if(is.null(smax_priors)==F){
      sm=sr_prior_mod(type='static',AC=TRUE,smax_priors=smax_priors)
    }
  }
  
  if(type=='rw'&par=='a'){
    if(is.null(smax_priors)==T){
      sm=sr_prior_mod(type='rw',par='a')
    }
    if(is.null(smax_priors)==F){
      sm=sr_prior_mod(type='rw',par='a',smax_priors=smax_priors)
    }
  }
  if(type=='rw'&par=='b'){
    if(is.null(smax_priors)==T){
      sm=sr_prior_mod(type='rw',par='b')
    }
    if(is.null(smax_priors)==F){
      sm=sr_prior_mod(type='rw',par='b',smax_priors=smax_priors)
    }
  }
  if(type=='rw'&par=='both'){
    if(is.null(smax_priors)==T){
      sm=sr_prior_mod(type='rw',par='both')
    }
    if(is.null(smax_priors)==F){
      sm=sr_prior_mod(type='rw',par='both',smax_priors=smax_priors)
    }
  }
  
  fit<-rstan::sampling(sm, data=datm,
                       control = control, warmup = warmup, 
                       chains = chains, iter = iter,verbose=FALSE)
  
  yrep=rstan::extract(fit,pars='y_rep')
  
  
  
}


#' Posterior predictive check plot
#'
#' @param fit A model fit in rstan from e.g. ricker_stan, ricker_rw_stan, ricker_hmm_stan
#' @param data the data that were used for the model in 'fit'
#' @export
post_check<- function(fit,data,pdf=FALSE,path=NULL,filename=NULL){
  if(is.null(fit$fit)==T){
    print('Model fit must specify: full_posterior=TRUE')
  }
  yrep_RS=rstan::extract(fit$fit,pars='y_rep',permuted=FALSE)
  yrep_R=array(NA,dim=dim(yrep_RS))
  for(t in 1:6){
    yrep_R[,t,]=log10(exp(yrep_RS[,t,])*data$S)
  }
  
  emp_RS=density(data$logRS,bw=0.02) #empirical distribution of log(R/S)
  emp_R=density(log10(data$R),bw=0.02) #empirical distribution of recruits
  
  cols=RColorBrewer::brewer.pal(n=7,'Blues')
  
  if(pdf==T){pdf(paste(path,'/',filename,'.pdf',sep=''),height=12,width=8)}
  par(mfrow=c(2,2))
  
  hist(data$logRS,main='',xaxt='n',breaks=30,freq=T,xlab='',ylab='counts of observations in time-series',col=adjustcolor('darkred',alpha.f=0.6),border='white',xlim=c(min(yrep_RS),max(yrep_RS)),xaxt='n')
  text('empirical obs.',x=par('usr')[2]-((par('usr')[2]-par('usr')[1])*0.2),y=par('usr')[4]-((par('usr')[4]-par('usr')[3])*0.1),col=adjustcolor('darkred',alpha.f=0.6))
  par(new=T)
  plot(rep(0,length(emp_RS$x))~emp_RS$x,type='n',col='darkred',ylab='',xlab='log(R/S)',xlim=c(min(yrep_RS),max(yrep_RS)),ylim=c(0,1),yaxt='n',bty='l')
   for(i in 1:6){
    d=density(yrep_RS[,i,],bw=0.05)
    lines(d$y/max(d$y)~d$x,col=cols[i+1])
    text(paste('pred. chain',i,sep=' '),x=par('usr')[2]-((par('usr')[2]-par('usr')[1])*0.2),y=par('usr')[4]-((par('usr')[4]-par('usr')[3])*(0.1+0.08*i)),col=cols[i+1])
  }
  hist(log10(data$R),main='',xaxt='n',breaks=30,freq=T,xlab='',ylab='counts of observations in time-series',col=adjustcolor('darkred',alpha.f=0.6),border='white',xlim=c(min(yrep_R),max(yrep_R)),xaxt='n')
  text('empirical obs.',x=par('usr')[2]-((par('usr')[2]-par('usr')[1])*0.2),y=par('usr')[4]-((par('usr')[4]-par('usr')[3])*0.1),col=adjustcolor('darkred',alpha.f=0.6))
  par(new=T)
  plot(rep(0,length(emp_R$x))~emp_R$x,type='n',col='darkred',ylab='',xlab='recuits (log10 axis)',xlim=c(min(yrep_R),max(yrep_R)),ylim=c(0,1),yaxt='n',xaxt='n',bty='l')
  for(i in 1:6){
    d=density(yrep_R[,i,],bw=0.02)
    lines(d$y/max(d$y)~d$x,col=cols[i+1])
    text(paste('pred. chain',i,sep=' '),x=par('usr')[2]-((par('usr')[2]-par('usr')[1])*0.2),y=par('usr')[4]-((par('usr')[4]-par('usr')[3])*(0.1+0.08*i)),col=cols[i+1])
  }
  pow <- c(min(round(log10(data$R)))-1):c(max(round(log10(data$R))))
  axis(1, col="black", at=seq(min(round(log10(data$R)))-1,max(round(log10(data$R))),by=1),tcl=-0.45, cex.axis=1.2,labels=10^pow)
  ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
  axis(1, log10(ticksat), col="black", labels=NA,
       tcl=-0.2, lwd=0, lwd.ticks=1)
  
  smaxs=rstan::extract(fit$fit,pars=c('prior_Smax','Smax'))
  hist(c(smaxs$prior_Smax/1e3),breaks=30,freq=T,xlim=c(0,max(c(smaxs[[1]]/1e3,smaxs[[2]]/1e3))),xlab='spawners (1000s)',col=adjustcolor('darkorange',alpha.f=0.5),border='white',main='')
  par(new=T)
  hist(c(smaxs$Smax/1e3),breaks=30,freq=T,xlim=c(0,max(c(smaxs[[1]]/1e3,smaxs[[2]]/1e3))),xlab='spawners (1000s)',col=adjustcolor('navy',alpha.f=0.5),border='white',main='',yaxt='n')
  text('Smax prior',x=par('usr')[2]-((par('usr')[2]-par('usr')[1])*0.2),y=par('usr')[4]-((par('usr')[4]-par('usr')[3])*0.1),col='darkorange')
  text('Smax posterior',x=par('usr')[2]-((par('usr')[2]-par('usr')[1])*0.2),y=par('usr')[4]-((par('usr')[4]-par('usr')[3])*0.2),col='navy')
  
  plot(c(data$R/1e3)~c(data$S/1e3),xlab='spawners (1000s)',ylab='recruits (1000s)',type='n',ylim=c(0,max(data$R/1e3)),xlim=c(0,max(data$S/1e3)*1.2),bty='l')
  sn=seq(0,max(data$S))
  muR=exp(median(fit$samples[,grepl('log_a',colnames(fit$samples))])-median(fit$samples[,grepl('b',colnames(fit$samples))])*sn)*sn
  lines(c(muR/1e3)~c(sn/1e3),lwd=3)
  dsmax=density(c(smaxs$Smax/1e3),bw=0.1)
  dsy=dsmax$y/max(dsmax$y)*par('usr')[4]
  lines(dsy~dsmax$x,col='navy',lwd=3)
  text(y=c(data$R/1e3),x=c(data$S/1e3),data$by,col='darkred',cex=0.8,font=2)
  text(x=max(data$S/1e3),y=max(data$R/1e3)*0.8,'Smax posterior',col='navy')
  if(pdf==T){
    dev.off()
    dev.off()
    }
}

sr_plot2=function(df,mod,title,make.pdf=FALSE,path,type=c('static','rw','hmm'),par=c('a','b','both'),form=c('stan','tmb'),ac=FALSE,sr_only=FALSE){
  if(type=='static'){ #static====
    x_new=seq(min(df$S),max(df$S),length.out=200)
    
    if(form=='stan'){
      post=as.data.frame(mod$draws(format='draws_matrix'))
      pred_df=data.frame(pred=exp(median(post$log_a)-median(post$b)*x_new)*x_new,x_new=x_new)
    }
    if(form=='tmb'){
      pred_df=data.frame(pred=exp(mod$alpha-mod$beta*x_new)*x_new,x_new=x_new)
    }
    plot=ggplot2::ggplot(df, aes(S, R)) +
      geom_line(data=pred_df,aes(x=x_new,y=pred),linewidth=1.3)+
      geom_point(aes(colour = by),size=2.5) +
      scale_colour_viridis_c(name='Year')+
      ggtitle(title)+
      xlab("Spawners") + 
      ylab("Recruits")+
      xlim(0, max(df$S))+
      ylim(0, max(df$R))+
      theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
            strip.text = element_text(face="bold", size=12),
            axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
  }
  if(type=='rw'){
    x_new=seq(min(df$S),max(df$S),length.out=200)
    by_q=round(quantile(df$by,seq(0,1,by=0.1)))
    if(par=='a'){ #rw alpha=====
      if(form=='stan'){
        log_a=as.data.frame(mod$draws(variables='log_a',format='draws_matrix'))
        b=as.data.frame(mod$draws(variables='b',format='draws_matrix'))
        pred_df=data.frame(x_new)
        for(n in 1:length(by_q)){
          pred_df[,1+n]=exp(median(log_a[,match(by_q[n],df$by)])-median(b$b)*x_new)*x_new
        }
        alpha_df=data.frame(by=seq(min(df$by),max(df$by)),med=apply(log_a,2,median),l90=apply(log_a,2,quantile,0.1),u90=apply(log_a,2,quantile,0.9))
        plot2=ggplot2::ggplot(alpha_df, aes(by,med)) +
          geom_line(aes(x=by,y=med),linewidth=1.3)+
          geom_point(aes(colour = by),size=4) +
          scale_colour_viridis_c(name='Year')+
          geom_ribbon(aes(ymin =l90, ymax =u90), alpha = 0.2)+
          ggtitle(title)+
          xlab("Year") + 
          ylab("log(Alpha)")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
        
      }
      
      if(form=='tmb'){
        pred_df=data.frame(x_new)
        for(n in 1:length(by_q)){
          pred_df[,1+n]=exp(mod$alpha[match(by_q[n],df$by)]-mod$beta*x_new)*x_new
        }
        alpha_df=data.frame(by=df$by,med=mod$alpha)
        
        plot2=ggplot2::ggplot(alpha_df, aes(by,med)) +
          geom_hline(yintercept=0,linetype='dashed')+
          geom_line(aes(x=by,y=med),linewidth=1.3)+
          geom_point(aes(colour = by),size=4) +
          scale_colour_viridis_c(name='Year')+
          ggtitle(paste(title))+
          xlab("Year") + 
          ylab("log(Alpha)")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
        
        
      }
      
      plot1=ggplot2::ggplot(df, aes(S, R)) +
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = by_q[1]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = by_q[2]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,4],colour = by_q[3]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,5],colour = by_q[4]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,6],colour = by_q[5]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,7],colour = by_q[6]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = by_q[1]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = by_q[2]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,4],colour = by_q[3]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,5],colour = by_q[4]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,6],colour = by_q[5]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,7],colour = by_q[6]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,8],colour = by_q[7]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,9],colour = by_q[8]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,10],colour = by_q[9]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,11],colour = by_q[10]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,12],colour = by_q[11]),linewidth=1.3)+
        geom_point(aes(colour = by),size=2.5) +
        scale_colour_viridis_c(name='Year')+
        ggtitle(title)+
        xlab("Spawners") + 
        ylab("Recruits")+
        xlim(0, max(df$S))+
        ylim(0, max(df$R))+
        theme_classic(14)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
      
      
      legend = cowplot::get_legend(plot1)
      
      plot_rw_a=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                                   plot2 + theme(legend.position="none"),
                                   ncol=2,nrow=1,labels=c("A","B"))
      plot=cowplot::plot_grid(plot_rw_a,legend,rel_widths = c(3,.25))
      if(sr_only==TRUE){plot=plot1}
    }
    if(par=='b'){ ###rw beta=====
      if(form=='stan'){
        post=as.data.frame(mod$draws(format='draws_matrix'))
        pred_df=data.frame(x_new)
        
        for(n in 1:length(by_q)){
          pred_df[,1+n]=exp(median(post$log_a)-median(post$b[,match(by_q[n],df$by)])*x_new)*x_new
        }
        
        beta_df=data.frame(by=df$by,med=apply(post$S_max,2,median),l90=apply(post$S_max,2,quantile,0.1),u90=apply(post$S_max,2,quantile,0.9))
        
        plot2=ggplot2::ggplot(beta_df, aes(by,med)) +
          geom_line(aes(x=by,y=med),linewidth=1.3)+
          geom_point(aes(colour = by),size=4) +
          scale_colour_viridis_c(name='Year')+
          geom_ribbon(aes(ymin =l90, ymax =u90), alpha = 0.2)+
          ggtitle(paste(title))+
          xlab("Year") + 
          ylab("Smax")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
        
      }
      if(form=='tmb'){
        pred_df=data.frame(x_new)
        for(n in 1:length(by_q)){
          pred_df[,1+n]=exp(mod$alpha-mod$beta[match(by_q[n],df$by)]*x_new)*x_new
        }
        beta_df=data.frame(by=df$by,med=mod$Smax)
        
        plot2=ggplot2::ggplot(beta_df, aes(by,med)) +
          geom_line(aes(x=by,y=med),linewidth=1.3)+
          geom_point(aes(colour = by),size=4) +
          scale_colour_viridis_c(name='Year')+
          ggtitle(paste(title))+
          xlab("Year") + 
          ylab("Smax")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold",size=14),axis.title = element_text(face="bold",size=14),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
        
        
      }
      
      
      plot1=ggplot2::ggplot(df, aes(S, R)) +
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = by_q[1]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = by_q[2]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,4],colour = by_q[3]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,5],colour = by_q[4]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,6],colour = by_q[5]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,7],colour = by_q[6]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = by_q[1]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = by_q[2]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,4],colour = by_q[3]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,5],colour = by_q[4]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,6],colour = by_q[5]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,7],colour = by_q[6]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,8],colour = by_q[7]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,9],colour = by_q[8]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,10],colour = by_q[9]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,11],colour = by_q[10]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,12],colour = by_q[11]),linewidth=1.3)+
        geom_point(aes(colour = by),size=2.5) +
        scale_colour_viridis_c(name='Year')+
        ggtitle(title)+
        xlab("Spawners") + 
        ylab("Recruits")+
        xlim(0, max(df$S))+
        ylim(0, max(df$R))+
        theme_classic(14)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold",size=14),axis.title = element_text(face="bold",size=14),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
      
      
      legend = cowplot::get_legend(plot1)
      
      plot_rw_b=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                                   plot2 + theme(legend.position="none"),
                                   ncol=2,nrow=1,labels=c("A","B"))
      plot=cowplot::plot_grid(plot_rw_b,legend,rel_widths = c(3,.25))
      if(sr_only==TRUE){plot=plot1}
    }
    if(par=='both'){ #rw alpha beta=====
      if(form=='stan'){
        post=as.data.frame(mod$draws(format='draws_matrix'))
        pred_df=data.frame(x_new)
        for(n in 1:length(by_q)){
          pred_df[,1+n]=exp(median(post$log_a[,match(by_q[n],df$by)])-median(post$b[,match(by_q[n],df$by)])*x_new)*x_new
        }
        
        alphabeta_df=data.frame(by=df$by,a_med=apply(post$log_a,2,median),a_l90=apply(post$log_a,2,quantile,0.15),a_u90=apply(post$log_a,2,quantile,0.85),b_med=apply(post$S_max,2,median),b_l90=apply(post$S_max,2,quantile,0.1),b_u90=apply(post$S_max,2,quantile,0.9))
        
        plot2=ggplot2::ggplot(alphabeta_df, aes(by,a_med)) +
          geom_line(aes(x=by,y=a_med),linewidth=1.3)+
          geom_point(aes(colour = by),size=3) +
          scale_colour_viridis_c(name='Year')+
          geom_ribbon(aes(ymin =a_l90, ymax =a_u90), alpha = 0.2)+
          ggtitle(paste(title))+
          xlab("Year") + 
          ylab("log(Alpha)")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
        
        plot3=ggplot2::ggplot(alphabeta_df, aes(by,b_med)) +
          geom_line(aes(x=by,y=b_med),linewidth=1.3)+
          geom_point(aes(colour = by),size=3) +
          scale_colour_viridis_c(name='Year')+
          geom_ribbon(aes(ymin =b_l90, ymax =b_u90), alpha = 0.2)+
          ggtitle(paste(title))+
          xlab("Year") + 
          ylab("Smax")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
        
      }
      if(form=='tmb'){
        pred_df=data.frame(x_new)
        for(n in 1:length(by_q)){
          pred_df[,1+n]=exp(mod$alpha[match(by_q[n],df$by)]-mod$beta[match(by_q[n],df$by)]*x_new)*x_new
        }
        alphabeta_df=data.frame(by=df$by,a_med=mod$alpha,b_med=mod$Smax)
        
        plot2=ggplot2::ggplot(alphabeta_df, aes(by,a_med)) +
          geom_line(aes(x=by,y=a_med),linewidth=1.3)+
          geom_point(aes(colour = by),size=3) +
          scale_colour_viridis_c(name='Year')+
          ggtitle(paste(title))+
          xlab("Year") + 
          ylab("log(Alpha)")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
        
        plot3=ggplot2::ggplot(alphabeta_df, aes(by,b_med)) +
          geom_line(aes(x=by,y=b_med),linewidth=1.3)+
          geom_point(aes(colour = by),size=3) +
          scale_colour_viridis_c(name='Year')+
          ggtitle(paste(title))+
          xlab("Year") + 
          ylab("Smax")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
        
        
      }
      
      plot1=ggplot2::ggplot(df, aes(S, R)) +
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = by_q[1]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = by_q[2]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,4],colour = by_q[3]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,5],colour = by_q[4]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,6],colour = by_q[5]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,7],colour = by_q[6]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = by_q[1]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = by_q[2]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,4],colour = by_q[3]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,5],colour = by_q[4]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,6],colour = by_q[5]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,7],colour = by_q[6]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,8],colour = by_q[7]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,9],colour = by_q[8]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,10],colour = by_q[9]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,11],colour = by_q[10]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,12],colour = by_q[11]),linewidth=1.3)+
        geom_point(aes(colour = by),size=2.5) +
        scale_colour_viridis_c(name='Year')+
        ggtitle(title)+
        xlab("Spawners") + 
        ylab("Recruits")+
        xlim(0, max(df$S))+
        ylim(0, max(df$R))+
        theme_classic(14)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
      
      
      legend = cowplot::get_legend(plot1)
      
      plot=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                              plot2 + theme(legend.position="none"),
                              plot3 + theme(legend.position="none"),
                              legend,
                              ncol=2,nrow=2,labels=c("A","B","C"))
      if(sr_only==TRUE){plot=plot1}
      
    }
    
  }
  if(type=='hmm'){ 
    
    x_new=seq(min(df$S),max(df$S),length.out=200)
    pred_df=data.frame(x_new)
    
    if(par=='a'){ #hmm alpha====
      
      if(form=='stan'){
        post=as.data.frame(mod$draws(format='draws_matrix'))
        pred_df[,2]=exp(median(post$log_a[,1])-median(post$b)*x_new)*x_new
        pred_df[,3]=exp(median(post$log_a[,2])-median(post$b)*x_new)*x_new
        df$gamma=apply(post$gamma[,,2],2,median)
        gamma_df=data.frame(by=df$by,gamma=apply(post$gamma[,,2],2,median),gamma_l90=apply(post$gamma[,,2],2,quantile,0.1),gamma_u90=apply(post$gamma[,,2],2,quantile,0.9))
        
        plot2=ggplot2::ggplot(gamma_df, aes(by,gamma)) +
          ylim(0,1)+
          geom_hline(yintercept=0.5,linetype='dashed')+
          geom_line(aes(x=by,y=gamma),linewidth=1.3)+
          geom_point(aes(colour = gamma),size=4) +
          scale_colour_viridis_c(name='p')+
          geom_ribbon(aes(ymin =gamma_l90, ymax =gamma_u90), alpha = 0.2)+
          ggtitle(paste(title))+
          xlab("Year") + 
          ylab("Prob. of high prod. regime")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
      }
      
      if(form=='tmb'){
        pred_df[,2]=exp(mod$alpha[1]-mod$beta*x_new)*x_new
        pred_df[,3]=exp(mod$alpha[2]-mod$beta*x_new)*x_new
        gamma_df=data.frame(by=df$by,gamma=mod$probregime[2,])
        
        plot2=ggplot2::ggplot(gamma_df, aes(by,gamma)) +
          ylim(0,1)+
          geom_hline(yintercept=0.5,linetype='dashed')+
          geom_line(aes(x=by,y=gamma),linewidth=1.3)+
          geom_point(aes(colour = gamma),size=4) +
          scale_colour_viridis_c(name='P(High prod. regime)')+
          ggtitle(paste(title))+
          xlab("Year") + 
          ylab("Prob. of high prod. regime")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
        
      }
      
      plot1=ggplot2::ggplot(df, aes(S, R)) +
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = min(gamma_df$gamma)),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = max(gamma_df$gamma)),linewidth=1.3)+
        geom_point(aes(colour = gamma_df$gamma),size=2.5) +
        scale_colour_viridis_c(name='p')+
        ggtitle(title)+
        xlab("Spawners") + 
        ylab("Recruits")+
        theme_classic(14)+
        xlim(0, max(df$S))+
        ylim(0, max(df$R))+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
      
      
      
      legend = cowplot::get_legend(plot1)
      
      plot_hmm_a=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                                    plot2 + theme(legend.position="none"),
                                    ncol=2,nrow=1,labels=c("A","B"))
      plot=cowplot::plot_grid(plot_hmm_a,legend,rel_widths = c(3,.3))
      if(sr_only==TRUE){plot=plot1}
    }
    
    if(par=='b'){ #hmm beta====
      
      if(form=='stan'){
        post=as.data.frame(mod$draws(format='draws_matrix'))
        pred_df[,2]=exp(median(post$log_a)-median(post$b[,1])*x_new)*x_new
        pred_df[,3]=exp(median(post$log_a)-median(post$b[,2])*x_new)*x_new
        df$gamma=apply(post$gamma[,,2],2,median)
        gamma_df=data.frame(by=df$by,gamma=apply(post$gamma[,,1],2,median),gamma_l90=apply(post$gamma[,,1],2,quantile,0.1),gamma_u90=apply(post$gamma[,,1],2,quantile,0.9))
        
        plot2=ggplot2::ggplot(gamma_df, aes(by,gamma)) +
          ylim(0,1)+
          geom_hline(yintercept=0.5,linetype='dashed')+
          geom_line(aes(x=by,y=gamma),linewidth=1.3)+
          geom_point(aes(colour = gamma),size=4) +
          scale_colour_viridis_c(name='p')+
          geom_ribbon(aes(ymin =gamma_l90, ymax =gamma_u90), alpha = 0.2)+
          ggtitle(paste(title))+
          xlab("Year") + 
          ylab("Prob. of high cap. regime")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
      }
      if(form=='tmb'){
        pred_df[,2]=exp(mod$alpha-mod$beta[1]*x_new)*x_new
        pred_df[,3]=exp(mod$alpha-mod$beta[2]*x_new)*x_new
        gamma_df=data.frame(by=df$by,gamma=mod$probregime[1,])
        
        plot2=ggplot2::ggplot(gamma_df, aes(by,gamma)) +
          ylim(0,1)+
          geom_hline(yintercept=0.5,linetype='dashed')+
          geom_line(aes(x=by,y=gamma),linewidth=1.3)+
          geom_point(aes(colour = gamma),size=4) +
          scale_colour_viridis_c(name='p')+
          ggtitle(paste(title))+
          xlab("Year") + 
          ylab("Prob. of high cap. regime")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
        
      }
      
      plot1=ggplot2::ggplot(df, aes(S, R)) +
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = max(gamma_df$gamma)),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = min(gamma_df$gamma)),linewidth=1.3)+
        geom_point(aes(colour = gamma_df$gamma),size=2.5) +
        scale_colour_viridis_c(name='p')+
        ggtitle(title)+
        xlab("Spawners") + 
        ylab("Recruits")+
        xlim(0, max(df$S))+
        ylim(0, max(df$R))+
        theme_classic(14)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
      
      
      
      legend = cowplot::get_legend(plot1)
      
      plot_hmm_b=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                                    plot2 + theme(legend.position="none"),
                                    ncol=2,nrow=1,labels=c("A","B"))
      plot=cowplot::plot_grid(plot_hmm_b,legend,rel_widths = c(3,.3))
      if(sr_only==TRUE){plot=plot1}
    }
    if(par=='both'){ #hmm alpha beta====
      if(form=='stan'){
        post=as.data.frame(mod$draws(format='draws_matrix'))
        pred_df[,2]=exp(median(post$log_a[,1])-median(post$b[,1])*x_new)*x_new
        pred_df[,3]=exp(median(post$log_a[,2])-median(post$b[,2])*x_new)*x_new
        df$gamma=apply(post$gamma[,,2],2,median)
        gamma_df=data.frame(by=df$by,gamma=apply(post$gamma[,,2],2,median),gamma_l90=apply(post$gamma[,,2],2,quantile,0.1),gamma_u90=apply(post$gamma[,,2],2,quantile,0.9))
        
        plot2=ggplot2::ggplot(gamma_df, aes(by,gamma)) +
          ylim(0,1)+
          geom_hline(yintercept=0.5,linetype='dashed')+
          geom_line(aes(x=by,y=gamma),linewidth=1.3)+
          geom_point(aes(colour = gamma),size=4) +
          scale_colour_viridis_c(name='p')+
          geom_ribbon(aes(ymin =gamma_l90, ymax =gamma_u90), alpha = 0.2)+
          ggtitle(paste(title))+
          xlab("Year") + 
          ylab("Prob. of high prod. regime")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
      }
      if(form=='tmb'){
        pred_df[,2]=exp(mod$alpha[1]-mod$beta[1]*x_new)*x_new
        pred_df[,3]=exp(mod$alpha[2]-mod$beta[2]*x_new)*x_new
        gamma_df=data.frame(by=df$by,gamma=mod$probregime[1,])
        
        plot2=ggplot2::ggplot(gamma_df, aes(by,gamma)) +
          ylim(0,1)+
          geom_hline(yintercept=0.5,linetype='dashed')+
          geom_line(aes(x=by,y=gamma),linewidth=1.3)+
          geom_point(aes(colour = gamma),size=4) +
          scale_colour_viridis_c(name='p')+
          ggtitle(paste(title))+
          xlab("Year") + 
          ylab("Prob. of high cap. regime")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
      }
      
      plot1=ggplot2::ggplot(df, aes(S, R)) +
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = min(gamma_df$gamma)),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = max(gamma_df$gamma)),linewidth=1.3)+
        geom_point(aes(colour = gamma_df$gamma),size=2.5) +
        scale_colour_viridis_c(name='p')+
        ggtitle(title)+
        xlab("Spawners") + 
        ylab("Recruits")+
        xlim(0, max(df$S))+
        ylim(0, max(df$R))+
        theme_classic(14)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
      
      legend = cowplot::get_legend(plot1)
      
      plot_hmm_ab=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                                     plot2 + theme(legend.position="none"),
                                     ncol=2,nrow=1,labels=c("A","B"))
      plot=cowplot::plot_grid(plot_hmm_ab,legend,rel_widths = c(3,.3))
      if(sr_only==TRUE){plot=plot1}
    }
  }
  if(make.pdf==TRUE){
    if(type=='static'&ac==FALSE){ pdf(here(path,paste(paste(title,type,form,sep='_'),'.pdf',sep='')),width=8,height=6)}
    if(type=='static'&ac==TRUE){ pdf(here(path,paste(paste(title,type,'ac',form,sep='_'),'.pdf',sep='')),width=8,height=6)}
    if(type=='rw'&par=='both'){ pdf(here(path,paste(paste(title,type,par,form,sep='_'),'.pdf',sep='')),width=10,height=10)}
    if(type=='rw'&par!='both'){ pdf(here(path,paste(paste(title,type,par,form,sep='_'),'.pdf',sep='')),width=14,height=6)}
    if(type=='hmm'&par=='both'){ pdf(here(path,paste(paste(title,type,par,form,sep='_'),'.pdf',sep='')),width=8,height=6)}
    if(type=='hmm'&par!='both'){ pdf(here(path,paste(paste(title,type,par,form,sep='_'),'.pdf',sep='')),width=14,height=6)}
    print(plot)
    dev.off()
  }
  if(make.pdf==FALSE){
    return(plot)
  }
  
}


#plot functions from Gottfried Pestal####

#' Predicted vs observed plot
#'
#' @param post.obj posterior samples
#' @param main.title plot title
#' @param plot.scale plot scale - how to rescale spawners (e.g. thousands, millions, etc.)
#' @param scale.label plot scale label
#' @param rec.label y-axis label - default 'Adult Recruits'
#' @param spn.label x-axis label - default 'Spawners'
#' @param file.name filename for output
#' @export

plot_predvsobs <- function(post.obj = NULL,
                          main.title = "Stock - Predicted vs. Observed",
                          #plot.log = FALSE, not yet implemented
                          plot.scale = 10^6,
                          scale.label = "(Millions)",
                          rec.label = "Adult Recruits",
                          spn.label = "Spawners",
                          file.name = "Plot_FittedvsObs.png"
){
  
  png(filename = file.name,
      width = 480*4, height = 480*4.5, 
      units = "px", pointsize = 14*3.3, bg = "white",  res = NA)
  
  layout(matrix(c(1,2),ncol=1),heights=c(1,1.8))
  par(mai=c(3.1,3.1,0.8,0.6))
  
  rec.max <- max(post.obj$quants.obs$R,post.obj$quants.obs$p90,na.rm=TRUE)/plot.scale
  spn.max <- max(post.obj$quants.obs$S,post.obj$quants.curve$S, na.rm=TRUE)/plot.scale
  
  if(post.obj$model_type == "TVP"){main.title <- paste0(main.title," (Last2 Gen)")}
  
  
  plot(post.obj$quants.obs$by,post.obj$quants.obs$R/plot.scale,ylim=c(0,rec.max),
       xlab="Brood Year",ylab=paste(rec.label,scale.label),las=1,
       pch=21,col="darkblue",type="p",bty="n", main = main.title )
  
  polygon(c(post.obj$quants.obs$by,rev(post.obj$quants.obs$by)),
          c(post.obj$quants.obs$p75/plot.scale, rev(post.obj$quants.obs$p25/plot.scale)  ),
          col="tomato",border="tomato")
  lines(post.obj$quants.obs$by,post.obj$quants.obs$p90/plot.scale,col="red",lty=2)
  lines(post.obj$quants.obs$by,post.obj$quants.obs$p10/plot.scale,col="red",lty=2)
  lines(post.obj$quants.obs$by,post.obj$quants.obs$p50/plot.scale,col="red",lwd=3)
  
  # obs points color fade
  # this will have different shading sequence for each stock depending on number of years
  n.obs <- length(min(post.obj$quants.obs$by):max(post.obj$quants.obs$by))
  #bg.fade <- c(rep(0,5),seq(0,1, length.out=n.obs-5 ))
  bg.fade <- seq(0.025,1, length.out=n.obs )
  bg.col <- rgb(0,0,1,bg.fade)
  points(post.obj$quants.obs$by,post.obj$quants.obs$R/plot.scale,pch=21,col="darkblue",
         bg = "white")
  points(post.obj$quants.obs$by,post.obj$quants.obs$R/plot.scale,pch=21,col="darkblue",
         bg = bg.col)
  
  
  
  
  plot(post.obj$quants.obs$S/plot.scale,post.obj$quants.obs$R/plot.scale,
       ylim=c(0,rec.max),xlim=c(0,spn.max),
       xlab=paste(spn.label,scale.label),ylab=paste(rec.label,scale.label),las=1,
       pch=21,col="darkblue",type="p",bty="n",cex=1.3)
  
  
  polygon(c(post.obj$quants.curve$Spn/plot.scale,rev(post.obj$quants.curve$Spn/plot.scale)),
          c(post.obj$quants.curve$p75/plot.scale, rev(post.obj$quants.curve$p25/plot.scale)  ),
          col="tomato",border="tomato")
  lines(post.obj$quants.curve$Spn/plot.scale,post.obj$quants.curve$p90/plot.scale,col="red",lty=2)
  lines(post.obj$quants.curve$Spn/plot.scale,post.obj$quants.curve$p10/plot.scale,col="red",lty=2)
  lines(post.obj$quants.curve$Spn/plot.scale,post.obj$quants.curve$p50/plot.scale,col="red",lwd=3)
  
  
  points(post.obj$quants.obs$S/plot.scale,post.obj$quants.obs$R/plot.scale,pch=21,col="darkblue",bg = "white",cex=1.3)
  points(post.obj$quants.obs$S/plot.scale,post.obj$quants.obs$R/plot.scale,pch=21,col="darkblue",
         bg = bg.col,cex=1.3)
  
  
  text(post.obj$quants.obs$S/plot.scale,post.obj$quants.obs$R/plot.scale,
       labels =  post.obj$quants.obs$by,col="darkgrey",cex=0.8,adj=-0.2, xpd=NA)
  
  
  dev.off()
  
  
  
  
} # end plot_predvsobs 


plot_posteriors <- function(post.obj = NULL,
                            main.title = "Stock - Posteriors",
                            plot.scale = 10^6,
                            scale.label = "(Millions)",
                            file.name = "Plot_Posteriors.png",
                            trim = 0  # proportion of distr to trim at upper tail for abd variables
                            #(e.g., trim = 0.05 excludes largest 5% of samples)
){
  
  png(filename = file.name,
      width = 480*4, height = 480*5, 
      units = "px", pointsize = 14*4, bg = "white",  res = NA)
  
  par(mfrow = c(3,2))
  par(mai=c(3.1,3.1,4,0.6))
  
  
  vals.plot <- post.obj$posteriors$log_a
  label.plot <- paste("log_a")
  hist(vals.plot,breaks=40,main=label.plot,
       axes=FALSE,xlab="",ylab="",border="lightblue",col="lightblue")
  abline(v= mean(vals.plot),col="red" );axis(1)
  
  vals.plot <- post.obj$posteriors$b
  label.plot <- paste("b")
  hist(vals.plot,breaks=40,main=label.plot,
       axes=FALSE,xlab="",ylab="",border="lightblue",col="lightblue")
  abline(v= mean(vals.plot),col="red" );axis(1)
  
  
  vals.plot <- post.obj$posteriors$S_gen
  label.plot <- paste("S_gen", scale.label)
  
  if(trim == 0){
    hist(vals.plot,breaks=40,main=label.plot,
         axes=FALSE,xlab="",ylab="",border="lightblue",col="lightblue")
    abline(v= mean(vals.plot),col="red" );axis(1)
  }
  
  if(trim > 0){
    
    vals.plot <- vals.plot[vals.plot < quantile(vals.plot,probs=1-trim,na.rm=TRUE)]
    
    hist(vals.plot,breaks=40,main=paste0(label.plot," trim=",trim*100,"%"),
         axes=FALSE,xlab="",ylab="",border="lightblue",col="lightblue")
    abline(v= mean(vals.plot),col="red" );axis(1)
  }
  
  
  
  vals.plot <- post.obj$posteriors$S_msy
  label.plot <- paste("S_msy", scale.label)
  if(trim == 0){
    hist(vals.plot,breaks=40,main=label.plot,
         axes=FALSE,xlab="",ylab="",border="lightblue",col="lightblue")
    abline(v= mean(vals.plot),col="red" );axis(1)
  }
  
  if(trim > 0){
    
    vals.plot <- vals.plot[vals.plot < quantile(vals.plot,probs=1-trim,na.rm=TRUE)]
    
    hist(vals.plot,breaks=40,main=paste0(label.plot," trim=",trim*100,"%"),
         axes=FALSE,xlab="",ylab="",border="lightblue",col="lightblue")
    abline(v= mean(vals.plot),col="red" );axis(1)
  }
  
  vals.plot <- post.obj$posteriors$S_max
  label.plot <- paste("S_max", scale.label)
  if(trim == 0){
    hist(vals.plot,breaks=40,main=label.plot,
         axes=FALSE,xlab="",ylab="",border="lightblue",col="lightblue")
    abline(v= mean(vals.plot),col="red" );axis(1)
  }
  
  if(trim > 0){
    
    vals.plot <- vals.plot[vals.plot < quantile(vals.plot,probs=1-trim,na.rm=TRUE)]
    
    hist(vals.plot,breaks=40,main=paste0(label.plot," trim=",trim*100,"%"),
         axes=FALSE,xlab="",ylab="",border="lightblue",col="lightblue")
    abline(v= mean(vals.plot),col="red" );axis(1)
  }
  
  vals.plot <- post.obj$posteriors$U_msy
  label.plot <- paste("U_msy (Rate)")
  hist(vals.plot,breaks=40,main=label.plot,
       axes=FALSE,xlab="",ylab="",border="lightblue",col="lightblue")
  abline(v= mean(vals.plot),col="red" );axis(1)
  
  if(post.obj$model_type == "TVP"){main.title <- paste0(main.title," (Last2 Gen)")}
  
  title(main=main.title,outer=TRUE,line=-2)
  
  dev.off()
  
}


plotJointDistr <- function(joint.in, x.label ="Var 1",y.label ="Var 2",x.lim = NULL,y.lim=NULL){
  
  med.x <- median(joint.in[[1]])
  med.y <- median(joint.in[[2]])
  
  # colors for the contour lines
  #library(RColorBrewer)
  #k <- 8
  #my.cols <- rev(brewer.pal(k, "RdYlBu"))
  
  par(pty="s")
  plot(joint.in[[1]],joint.in[[2]], bty="n", xlab = x.label,ylab = y.label, las = 1,
       pch =19,col = "lightgrey",cex=0.7,cex.axis=1.4,cex.lab = 1.4,xlim=x.lim,ylim=y.lim)
  
  
  z <- kde2d(joint.in[[1]],joint.in[[2]], n=50)
  
  #contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE,lwd=4)
  contour(z, drawlabels=FALSE, nlevels=8, col="darkgrey", add=TRUE,lwd=4)
  
  
  
  x.dist <- par("usr")[2]-par("usr")[1]
  y.dist <- par("usr")[4]-par("usr")[3]
  
  RapidRicker::kernel_margin(joint.in[[2]],
                             at = par("usr")[2] - x.dist*0.15 ,width = x.dist*0.20,
                             col = "darkgrey",dir = "v")  
  
  RapidRicker::kernel_margin(joint.in[[1]],
                             at = par("usr")[4] - y.dist*0.13 ,width = y.dist*0.20,
                             col = "darkgrey",dir = "h")  
  
  abline(v = med.x,col="darkblue",lty=2,lwd=2)
  abline(h = med.y,col="darkblue",lty=2,lwd=2)
  
  
  
  RapidRicker::box_add(quantile(joint.in[[2]],probs=c(0.1,0.25,0.5,0.75,0.9)),
                       at = par("usr")[2] - x.dist*0.15 ,width = x.dist*0.05,
                       col = "darkblue",bg="lightgray",label="",dir = "v")
  
  RapidRicker::box_add(quantile(joint.in[[1]],probs=c(0.1,0.25,0.5,0.75,0.9)),
                       at = par("usr")[4] - y.dist*0.13 ,width = y.dist*0.05,
                       col = "darkblue",bg="lightgray",label="",dir = "h")
  
} #end plotJointDistr()


plot_jointpost <- function(joint.in = NULL,
                           x.label ="Var 1",y.label ="Var 2",x.lim = NULL,y.lim=NULL,
                           main.title = "Test Stock - Joint Posterior",
                           file.name = "Plot_FittedvsObs.png"
                           
){
  
  png(filename = file.name,
      width = 480*4, height = 480*4, 
      units = "px", pointsize = 14*3.3, bg = "white",  res = NA)
  
  par(mai=c(3.3,3.3,4,0.6))
  
  plotJointDistr(joint.in = joint.in, 
                 x.label = x.label,
                 y.label = y.label,
                 x.lim = x.lim,y.lim=y.lim)
  
  title(main=main.title)
  
  dev.off()		   
  
}		


plot_spnvsbm <- function(post.obj = NULL,
                         spn.data.all = NULL,
                         main.title = "Test Stock - Relative Abundance Benchmarks",
                         plot.scale = 10^6,
                         scale.label = "(Millions)",
                         AvgGen = 4 , # yrs to use for running avg (feed in dominant age class for each stock)
                         spn.label = "Spawners",
                         file.name = "Plot_RelAbdBM.png"
){
  
  png(filename = file.name,
      width = 480*4.5, height = 480*3.7, 
      units = "px", pointsize = 14*3.5, bg = "white",  res = NA)
  
  #layout(matrix(c(1,2),ncol=1),heights=c(1,1.8))
  #par(mai=c(3.1,3.1,0.8,0.6))
  
  if(post.obj$model_type == "TVP"){main.title <- paste0(main.title," (Last2 Gen)")}
  
  
  bm.df <- post.obj$quants.post %>% dplyr::filter(Variable %in% c("S_gen","S_msy","S_max")) %>%
    dplyr::select(Variable,p10,p25,p50,p75,p90)
  
  sr.used <- post.obj$source_data
  
  spn.in <- spn.data.all
  
  # fill in missing yrs with NA(to get proper line gaps, and for running geomean calc)
  full.yrs <-  min(spn.in$Year):max(spn.in$Year)
  missing.yrs <-  setdiff(full.yrs,spn.in$Year)
  #print(missing.yrs)
  if(length(missing.yrs>0)){ spn.in <- left_join(data.frame(by=full.yrs),spn.in,by="Year") }
  
  ylim <- c(0, max(spn.in$Spn,na.rm=TRUE))  # ADD BM RANGES!!!!!!!!!!
  #print(ylim)
  
  axis.scale <- 1
  axis.scale.label <- ""
  
  if(ylim[2] >= 1000 & ylim[2] < 10000){ axis.scale <- 100; axis.scale.label <- "(100s)"   }
  if(ylim[2] >= 10000 & ylim[2] < 1000000){ axis.scale <- 1000; axis.scale.label <- "(1000s)"   }
  if(ylim[2] >= 1000000 ){ ylim[2]<- 1000000; axis.scale.label <- "(Mill)"   }
  
  spn.in$Spn <- spn.in$Spn/axis.scale
  
  y.label.use <- paste("Spn",axis.scale.label)
  
  x.range <- range(spn.in$Year,na.rm=TRUE)
  x.ticks <- pretty(x.range)
  x.ticks <- x.ticks[x.ticks >=min(x.range) & x.ticks <= max(x.range)]
  
  
  plot(1:5,1:5,type="n",bty="n",ylim = ylim/axis.scale ,xlim = x.range+c(0,30) ,  las=1,
       xlab="", ylab = y.label.use,axes=FALSE, main=main.title, col.main = "darkblue")
  axis(2,las=1)
  axis(1,x.ticks)
  
  # calculate running geomean
  gm.in <- log(spn.in$Spn)
  
  
  gm.out <- exp(stats::filter(gm.in,rep(1/AvgGen,times =AvgGen),sides = 1))# 
  gm.out
  
  
  
  sgen.vec <- bm.df %>% dplyr::filter(Variable == "S_gen") %>% select(-Variable) %>% unlist()
  smax.vec <- bm.df %>% dplyr::filter(Variable == "S_max") %>% select(-Variable) %>% unlist()
  smsy.vec <- bm.df %>% dplyr::filter(Variable == "S_msy") %>% select(-Variable) %>% unlist()
  #print(sgen.vec)
  
  abline(h=sgen.vec[3]/axis.scale,col="red",lwd=3)
  abline(h=smsy.vec[3]*0.8/axis.scale,col="green",lwd=3,lty=2)
  
  # LOWER BM
  
  offset.use <- 3
  box_add(sgen.vec/axis.scale,at=x.range[2]+offset.use,width=2,col="darkblue",bg="lightblue",label="")
  box_add(smax.vec*0.2/axis.scale,at=x.range[2]+offset.use+4,width=2,col="darkblue",bg="lightblue",label="")
  axis(1,at=x.range[2]+c(offset.use,offset.use+4),labels=c("Sgen","20%Smax"),las=2)
  
  
  # UPPER bm
  
  offset.use <- 13
  box_add(smsy.vec*0.8/axis.scale,at=x.range[2]+offset.use,width=2,col="darkblue",bg="lightblue",label="")
  box_add(smax.vec*0.4/axis.scale,at=x.range[2]+offset.use+4,width=2,col="darkblue",bg="lightblue",label="")
  axis(1,at=x.range[2]+c(offset.use,offset.use+4),labels=c("80%Smsy","40%Smax"),las=2)
  
  # BIOL BM
  offset.use <- 23
  box_add(smsy.vec/axis.scale,at=x.range[2]+offset.use,width=2,col="darkblue",bg="lightblue",label="")
  box_add(smax.vec/axis.scale,at=x.range[2]+offset.use+4,width=2,col="darkblue",bg="lightblue",label="")
  axis(1,at=x.range[2]+c(offset.use,offset.use+4),labels=c("Smsy","Smax"),las=2)
  
  
  lines(spn.in$Year,spn.in$Spn,type="o",col="darkblue",pch=21,bg="white",cex=0.8)
  points(sr.used$by,sr.used$S/axis.scale,col="darkblue",pch=19,cex=0.8)
  
  lines(spn.in$Year,gm.out,type="l",col="red",lwd=4)
  
  
  legend("topleft",legend = c("Used for SR Fit", "Not Used", paste0(AvgGen,"yr Running GeoMean")),
         pch = c(19,21,NA), col= c("darkblue","darkblue","red"), lty=c(NA,NA,1), lwd=c(NA,NA,3),cex=0.7, bty="n" )
  
  
  
  dev.off()
  
  
  
  
} # end plot_spnvsbm



