function step_mail(n,EarthFaultTime)
setpref('Internet','E_mail','yuyo.jushin@gmail.com');
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username','yuyo.jushin@gmail.com');
setpref('Internet','SMTP_Password','70158230');
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
%str1 = num2str(loopnum);
str = ['macmini EFT ',num2str(EarthFaultTime),' n ',num2str(n)];
%strcat(info,str1);
sendmail('step-by-step.t-a@ezweb.ne.jp','Matlab',str)