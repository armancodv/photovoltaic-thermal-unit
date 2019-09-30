program pvt
    implicit none
    integer,parameter :: imax=101,jmax=86
    double precision,parameter :: pi=3.141592654,convergency_i=10.0**-6.0,convergency_t=10.0**-9.0
    double precision :: h,l,ratio,dt,final_time,time,sigma,v_inf,tamb,tsky,tsun,rho_g,k_g,c_g,alpha_g,eps_g,tau_g,d_g,k_a,v_a,beta_a,d_a,pac,beta,eta0,eps_pv,rho_pv,k_pv,c_pv,alpha_pv,d_pv,rho_ab,c_ab,d_ab,k_ab,diameter,per,a_t,a_t_ab,a_t_in,rho_t,c_t,k_t,d_t,rho_f,rho_n,rho_nf,c_f,c_n,c_nf,m_f,a_f,tf_inlet,nu_f,k_f,k_n,k_nf,q0,const,hconvt_f,rho_i,c_i,k_i,d_i,ra,nu_a,zz1,zz2,delta_i_g,delta_i_pv,delta_i_ab,delta_i_i,delta_i_t,delta_i_f,delta_i,delta_t_g,delta_t_pv,delta_t_ab,delta_t_i,delta_t_t,delta_t_f,delta_t,k_ad,d_ad,h_g_pv,h_pv_ab,h_ab_t,h_ab_i,h_t_i,k1,k2,x_pipe,x_center,y_pipe,x_margin,y_margin,total_residual,dx_t,dy_t,g,h_wind,a,b,c,dw,de,dn,ds,dv_t,dv_f,gamma,phi,etha_pvt,etha_th,etha_el,etha_tht,e_th,e_in,e_el,ex_th,ex_in,ex_el,eps_pvt,eps_th,eps_el,tave,d_n,begin_time,end_time,l_pos,mu_nf,pr_nf,re_nf,f_f,mu_f,v_f,k_beta,k_a1,k_a2,k_a3,en_f_t_total,en_f_t2_total,en_f_f_total,f_nf,phiw,en_system,k_fx,dp,dp_dv
    double precision :: h_pv_te, h_te_ab, e_te, etha_el_te, etha_el_pv, e_el_pv, e_el_te
    double precision :: p_s, p_s_g, p_s_pv, p_g_ref, p_pv, p_fluid, p_g_conv, p_g_pv, p_g_sky, p_pv_ab, p_ab_t, p_ab_in, p_t_in, p_in_conv
    integer :: i,j,k,iterate,counter,time_iter,n_pipe,i_margin,j_margin,i_pipe,i_center,j_pipe,i_pos,j_pos,pipe_counter,element_counter,i_radius,j_radius,total_element,i_marginf,j_marginf,file,file_length,i_length,i_file
    double precision,dimension(5000) :: l_line
    double precision,dimension(imax) :: x,dx_pe,dx_wp,dx
    double precision,dimension(jmax) :: y,dy_pn,dy_sp,dy
    double precision,dimension(imax,jmax) :: tg,tpv,tab,tt,tf,ti,tg_old,tpv_old,tab_old,ti_old,tg_pre,tpv_pre,tab_pre,ti_pre,sp,sc,ae,aw,an,as,ap,ap0,ke,kw,kn,ks,pipe_line,pipe_shape,dv,en_f_t,en_f_f
    double precision,dimension(imax,jmax) :: tte
    double precision,dimension(:),allocatable :: kr,kl,al,ar,apt,ap0t,spt,sct,tt1,tf1,tt1_old,tf1_old,tt1_pre,tf1_pre,apf,ap0f,spf,scf,fr,fl,dr,dl,pr,pl,en1_t,en1_f
    double precision,dimension(:),allocatable :: phiw_array,g_array,tamb_array,tf_array,ex_outlet,ex_surface,ex_elef,test_time,m_f_array
    character(11),dimension(:),allocatable :: folder_array
    character(8) :: date
    character(4) :: prefix
    character(5) :: prefix2
    character(2) :: folder

    ! Al2O3 0.002wt%
    i_length=13
    file_length=i_length*10
    allocate(folder_array(file_length),phiw_array(file_length),g_array(file_length),tamb_array(file_length),tf_array(file_length),ex_outlet(file_length),ex_surface(file_length),ex_elef(file_length),test_time(file_length),m_f_array(file_length))
    rho_n=3970.0     !kg/m^3
    c_n=765.0         !j/kgk (nano)
    k_n=40.0           !w/mK
    d_n=2.0*10.0**(-8.0)
    prefix='PvAl'
    phiw_array=0.002
    ex_outlet=0.0
    ex_surface=0.0
    ex_elef=0.0
    do file=1,file_length
        i_file = mod(file,i_length)
        if (i_file==0) then
            i_file = i_length
        end if
        test_time(file)=0.0
        if (int((file-1)/i_length+1)==1) then
            g_array(file)=620.0
            tf_array(file)=307.0+(316.0-307.0)*real(i_file-1)/real(i_length-1)
            m_f_array(file)=0.00834
            tamb_array(file)=306.0
            prefix2='TfMin'
        else if (int((file-1)/i_length+1)==2) then
            g_array(file)=920.0
            tf_array(file)=307.0+(316.0-307.0)*real(i_file-1)/real(i_length-1)
            m_f_array(file)=0.00834
            tamb_array(file)=306.0
            prefix2='TfAve'
        else if (int((file-1)/i_length+1)==3) then
            g_array(file)=1050.0
            tf_array(file)=307.0+(316.0-307.0)*real(i_file-1)/real(i_length-1)
            m_f_array(file)=0.00834
            tamb_array(file)=306.0
            prefix2='TfMax'
        else if (int((file-1)/i_length+1)==4) then
            g_array(file)=620.0+(1050.0-620.0)*real(i_file-1)/real(i_length-1)
            tf_array(file)=307.0
            m_f_array(file)=0.00834
            tamb_array(file)=306.0
            prefix2='SoMin'
        else if (int((file-1)/i_length+1)==5) then
            g_array(file)=620.0+(1050.0-620.0)*real(i_file-1)/real(i_length-1)
            tf_array(file)=313.0
            m_f_array(file)=0.00834
            tamb_array(file)=306.0
            prefix2='SoAve'
        else if (int((file-1)/i_length+1)==6) then
            g_array(file)=620.0+(1050.0-620.0)*real(i_file-1)/real(i_length-1)
            tf_array(file)=316.0
            m_f_array(file)=0.00834
            tamb_array(file)=306.0
            prefix2='SoMax'
        else if (int((file-1)/i_length+1)==7) then
            g_array(file)=920.0
            tf_array(file)=313.0
            m_f_array(file)=0.001+(0.01-0.001)*real(i_file-1)/real(i_length-1)
            tamb_array(file)=306.0
            prefix2='MfLam'
        else if (int((file-1)/i_length+1)==8) then
            g_array(file)=920.0
            tf_array(file)=313.0
            m_f_array(file)=0.014+(0.02-0.014)*real(i_file-1)/real(i_length-1)
            tamb_array(file)=306.0
            prefix2='MfTur'
        else if (int((file-1)/i_length+1)==9) then
            g_array(file)=920.0
            tf_array(file)=313.0
            m_f_array(file)=0.00834
            tamb_array(file)=303.0+(309.0-303.0)*real(i_file-1)/real(i_length-1)
            prefix2='TiAve'
        else if (int((file-1)/i_length+1)==10) then
            g_array(file)=-9.3432*i_file**2.0+146.16*i_file+484.34
            tf_array(file)=-0.1314*i_file**2.0+2.4699*i_file+304.68
            m_f_array(file)=0.00834
            tamb_array(file)=0.5*i_file+302.5
            test_time(file)=0.5*i_file+9.0
            prefix2='DayEx'
        end if
        
        write(folder,'(i2)') i_file
        if(i_file<10) then
            folder(:1)='0';
        end if
        folder_array(file) = prefix//prefix2//folder
    end do
    
    

    call input_data
    call mesh_generation
    g=g_array(1)
    phiw=phiw_array(1)
    tamb=tamb_array(1)
    tf_inlet=tf_array(1)
    call initialize
    open(200,file=prefix//'_results.csv')
    open(201,file=prefix//'_plot_phi.plt')
    open(202,file=prefix//'_min.csv')
    open(203,file=prefix//'_powers.csv')
    write(200,'(a1723)') 'No,Test Time,Solar Irradiation,Ambient Temperature,Inlet Temperature,Nanofluid Mass Fraction,Nanofluid Volume Fraction,,Outlet Temperature,Glass Mean Temperature,PV Mean Temperature,TE Mean Temperature,Absorber Mean Temperature,Tube Mean Temperature,Fluid Mean Temperature,Insulation Mean Temperature,TE Temperature Difference,,Glass Max Temperature,PV Max Temperature,TE Max Temperature,Absorber Max Temperature,Tube Max Temperature,Fluid Max Temperature,Insulation Max Temperature,,Experiment Outlet Temperature,Experiment Surface Temperature,Diffrence Outlet Temperature,Diffrence Surface Temperature,,Overal Efficiency,Thermal Efficiency,PV Electrical Efficiency,TE Electrical Efficiency,Electrical Efficiency,Experiment Electrical Efficiency,Electrical Efficiency Diffrence,Overal Equivalent Efficiency,,Thermal Exergetic Efficiency,Electrical Exergetic Efficiency,Exergy Thermal Power,Exergy Solar irradiation power,Exergy Output electrical power,Entropy Generation (Thermal),Entropy Generation (Thermal-Total),Pressure Diffrence,Entropy Generation (Viscous),Entropy Generation (System),,Reynolds Number,Prandtl Number,Nusselt Number,,Heat Transfer Coefficient,Tube-Insulation Heat Transfer Coefficient,Absorber-Tube Heat Transfer Coefficient,Absorber-Insulation Heat Transfer Coefficient,PV-TE Heat Transfer Coefficient,TE-Absorber Heat Transfer Coefficient,Glass-PV Heat Transfer Coefficient,,Solar,Sun->Glass,Glass->PV,Sun->PV,PV Electrical,PV->TE,TE electrical,TE->Abs,Abs->Tube,Fluid Thermal,,Nanofluid Density,Nanofluid Specific Heat Capacity,Nanofluid Thermal conductivity,Nanofluid Viscosity,,Sky Temperature,Velocity,Mass Flow,,Date,During Time,Time Step,Final Time,X-direction Elements,Y-direction Elements'
    write(200,'(a249)') ',h,w/m^2k,K,K,%,%,,K,K,K,K,K,K,K,K,K,,K,K,K,K,K,K,K,,K,K,K,K,,%,%,%,%,%,%,%,%,,%,%,w/m^2,w/m^2,w/m^2,w/k,w/k,p,w/k,w/k,,,,,,W/(m2.K),W/(m2.K),W/(m2.K),W/(m2.K),W/(m2.K),W/(m2.K),W/(m2.K),,W,W,W,W,W,W,W,W,W,W,,kg/m^3,j/kg,w/mk,,,K,m/s,kg/s,,,s,s,s,,'
    write(200,'(a1054)') 'folder_array(file),test_time(file),g,tamb,tf_inlet,phiw*100.0,phi*100.0,,tf1(total_element),sum(tg)/size(tg),sum(tpv)/size(tpv),sum(tte)/size(tte),sum(tab)/size(tab),sum(tt1)/size(tt1),sum(tf1)/size(tf1),sum(ti)/size(ti),sum(tpv)/size(tpv)-sum(tab)/size(tab),,maxval(tg),maxval(tpv),maxval(tte),maxval(tab),maxval(tt1),maxval(tf1),maxval(ti),,ex_outlet(file),ex_surface(file),tf1(total_element)-ex_outlet(file),ex_surface(file)-maxval(tpv),,etha_pvt,etha_th,etha_el_pv,etha_el_te,etha_el,ex_elef(file)*100.0,ex_elef(file)*100.0-etha_el,etha_tht,,eps_th,eps_el,ex_th,ex_in,ex_el,en_f_t_total,en_f_t2_total,dp,en_f_f_total,en_system,,re_nf,pr_nf,nu_f,,hconvt_f,h_t_i,h_ab_t,h_ab_i,h_pv_te,h_te_ab,h_g_pv,,g*l*h,alpha_g*g*l*h,h_g_pv*(sum(tg)/size(tg)-sum(tpv)/size(tpv))*l*h,g*tau_g*l*h,e_el_pv*l*h,h_pv_te*(sum(tpv)/size(tpv)-sum(tte)/size(tte))*l*h,e_el_te*l*h,h_te_ab*(sum(tte)/size(tte)-sum(tab)/size(tab))*l*h,h_ab_t*(sum(tab)/size(tab)-sum(tt1)/size(tt1))*l*h,e_th*l*h,,rho_nf,c_nf,k_nf,mu_nf,,tsky,v_f,m_f,,date,end_time-begin_time,dt,time,imax,jmax'
    write(201,'(a650)')'variables="Volume Fraction (%)","Mass Fraction (%)","Test Time (h)","Outlet Temperature","PV Temperature","Experiment Outlet Temperature","Experiment PV Temperature","Overal Efficiency","Thermal Efficiency","Electrical Efficiency","Experiment Electrical Efficiency","Overal Equivalent Efficiency","Overal Exergetic Efficiency","Thermal Exergetic Efficiency","Electrical Exergetic Efficiency","Nanofluid Density","Nanofluid Specific Heat Capacity","Nanofluid Thermal conductivity","Nanofluid Viscosity","Reynolds Number","Prandtl Number","Nusselt Number","Entropy Generation (Thermal)","Entropy Generation (Thermal-Total)","Entropy Generation (Viscous)","Entropy Generation (System)"'
    write(202,'(a276)') 'Solar Irradiation,Inlet Temperature,Mass Flow Rate,Outlet Temperature,TE Temperature Difference,Overal Efficiency,Thermal Efficiency,PV Electrical Efficiency,TE Electrical Efficiency,Electrical Efficiency,Thermal Exergetic Efficiency,PV Electrical,PV->TE,TE electrical,Thermal'
    write(203,'(a114)') 'p_s, p_s_g, p_s_pv, p_g_ref, p_pv, p_fluid, p_g_conv, p_g_pv, p_g_sky, p_pv_ab, p_ab_t, p_ab_in, p_t_in, p_in_conv'

    do file=79,104
        write(*,*) '-------------------------------------------------'
        write(*,*) 'File=',file
        m_f=m_f_array(file)
        g=g_array(file)                  !w/m^2k
        phiw=phiw_array(file)
        tamb=tamb_array(file)
        tf_inlet=tf_array(file)

        call open_files
        call input_property

        call cpu_time(begin_time)
        iterate=0
        time_iter=0
    
        ! time loop
        do while (time<500000)
            time=time+dt
            time_iter=time_iter+1
            do counter=1,10000
                iterate=iterate+1
                call fluid_property
    
                ! glass
                call coefficients_glass
                do j=1,jmax
                    call tdma(aw(1:imax,j),ap(1:imax,j),ae(1:imax,j), sc(1:imax,j),imax,tg(1:imax,j))
                end do
    
                ! pv
                call coefficients_pv
                do j=1,jmax
                    call tdma(aw(1:imax,j),ap(1:imax,j),ae(1:imax,j), sc(1:imax,j),imax,tpv(1:imax,j))
                end do
    
                ! absorber
    	        call coefficients_absorber
	            do j=1,jmax
                    call tdma(aw(1:imax,j),ap(1:imax,j),ae(1:imax,j), sc(1:imax,j),imax,tab(1:imax,j))
                end do

                ! tube
                call coefficients_tube
                call tdma(al,apt,ar, sct,total_element,tt1)
	            call res_tube
     
                ! fluid
	            call coefficients_fluid
                call tdma(al,apf,ar,scf,total_element,tf1)
	            call res_fluid
    	 
                ! insulation
                call coefficients_insulation
	            do j=1,jmax
                    call tdma(aw(1:imax,j),ap(1:imax,j),ae(1:imax,j),sc(1:imax,j),imax,ti(1:imax,j))
                end do
            
                delta_i_g=maxval(abs((tg-tg_pre)/tg))
                delta_i_pv=maxval(abs((tpv-tpv_pre)/tpv))
                delta_i_ab=maxval(abs((tab-tab_pre)/tab))
                delta_i_t=maxval(abs((tt1-tt1_pre)/tt1))
                delta_i_f=maxval(abs((tf1-tf1_pre)/tf1))
                delta_i_i=maxval(abs((ti-ti_pre)/ti))
                delta_i=max(delta_i_g,delta_i_pv,delta_i_ab,delta_i_t,delta_i_f,delta_i_i)
                write(*,'(a11,a13,i6,a8,i6,a9,i6,a9,f14.12)') folder_array(file),' | Time Iter=',time_iter,' | Time=',nint(time),'s | Iter=',iterate,' | Error=',delta_i
                write(1,*) iterate,delta_i
                tg_pre=tg
                tpv_pre=tpv
                tab_pre=tab
                tt1_pre=tt1
                tf1_pre=tf1
                ti_pre=ti
                if(delta_i<=convergency_i) exit
            end do
            delta_t_g=maxval(abs((tg-tg_old)/tg))
            delta_t_pv=maxval(abs((tpv-tpv_old)/tpv))
            delta_t_ab=maxval(abs((tab-tab_old)/tab))
            delta_t_t=maxval(abs((tt1-tt1_old)/tt1))
            delta_t_f=maxval(abs((tf1-tf1_old)/tf1))
            delta_t_i=maxval(abs((ti-ti_old)/ti))
            delta_t=max(delta_t_g,delta_t_pv,delta_t_ab,delta_t_t,delta_t_f,delta_t_i)
            write(2,*) time_iter,delta_t
            tg_old=tg
            tpv_old=tpv
            tab_old=tab
            tt1_old=tt1
            tf1_old=tf1
            ti_old=ti
            if(delta_t<=convergency_t) exit
        end do
        call efficiency
        call entropy_fluid
        call energy_balance
        call output
    end do
    pause
contains

subroutine input_data
    l=0.63             !m along x axis==>x(i)
    h=0.54             !m along y axis==>y(j)
    diameter=0.008     !m diameter of cooling pipeline network
    n_pipe=16          !number of pipes
    x_pipe=0.538       !m length of pipe along l
    y_pipe=0.48       !m length of pipe along h
    x_center=0.382
    dt=20.0             !second
    dx_t=l/(imax-2.0)
    dy_t=h/(jmax-2.0)
    x_margin=(l-x_pipe)/2.0              !m margin along x axis or l=length
    y_margin=(h-y_pipe)/2.0       !m margin along y axis or h=height
    return
end subroutine

subroutine mesh_generation
    do i=1,imax
        x(i)=x_func(i)
    end do
    
    do j=1,jmax
        y(j)=y_func(j)
    end do

    do i=2,imax-1
        dx_pe(i)=dx_pe_func(i)
        dx_wp(i)=dx_wp_func(i)
        dx(i)=dx_func(i)
    end do
    dx_pe(1)=dx_pe_func(1)
    dx_wp(imax)=dx_wp_func(imax)
    dx(1)=dx_func(1)
    dx(imax)=dx_func(imax)

    do j=2,jmax-1
        dy_pn(j)=dy_pn_func(j)
        dy_sp(j)=dy_sp_func(j)
        dy(j)=dy_func(j)
    end do
    dy_pn(1)=dy_pn_func(1)
    dy_sp(jmax)=dy_sp_func(jmax)
    dy(1)=dy_func(1)
    dy(jmax)=dy_func(jmax)
    
    do i=2,imax-1
        do j=2,jmax-1
            dv(i,j)=dv_func(i,j)
        end do
    end do

    ! locate pipe line
    pipe_shape=0.0
    pipe_line=0.0

    i_margin=nint((x_margin-dx_t/2.0)/dx_t)+1 ! first x-element number of pipe 
    j_margin=nint((y_margin-dy_t/2.0)/dy_t)+1 ! first y-element number of pipe
    i_pipe=nint(x_pipe/dx_t) ! number of elements of pipe x-lendgh
    i_center=nint(x_center/dx_t) ! number of elements of pipe x-lendgh
    j_pipe=nint(y_pipe/(n_pipe*dy_t)) ! number of elements of pipe y-lendgh
    i_radius=nint((diameter/2.0)/dx_t) ! number of elements of x-diameter
    j_radius=nint((diameter/2.0)/dy_t) ! number of elements of y-diameter
    
    i_pos=i_margin ! x-element number
    j_pos=j_margin ! y-element number
    pipe_counter=1 ! pipe counter
    element_counter=0 ! element counter
    l_pos=-dx_t/2.0
    do
        pipe_counter=pipe_counter+1
        if(pipe_counter==9) then
            do i=i_pos,i_pos-i_center,-1
	            element_counter=element_counter+1
                pipe_line(i,j_pos)=1.0
	            pipe_shape(i,j_pos:j_pos+j_radius)=1.0
	            pipe_shape(i,j_pos-j_radius:j_pos)=1.0
                l_pos=l_pos+dx_t
                l_line(element_counter)=l_pos
            end do
        else if(pipe_counter==10) then
            do i=i_pos,i_pos+i_center,1
	            element_counter=element_counter+1
                pipe_line(i,j_pos)=1.0
	            pipe_shape(i,j_pos:j_pos+j_radius)=1.0
	            pipe_shape(i,j_pos-j_radius:j_pos)=1.0
                l_pos=l_pos+dx_t
                l_line(element_counter)=l_pos
            end do
        else if(mod(pipe_counter,2)==0) then
            do i=i_pos,i_pos+i_pipe,1
	            element_counter=element_counter+1
                pipe_line(i,j_pos)=1.0
	            pipe_shape(i,j_pos:j_pos+j_radius)=1.0
	            pipe_shape(i,j_pos-j_radius:j_pos)=1.0
                l_pos=l_pos+dx_t
                l_line(element_counter)=l_pos
            end do
        else
            do i=i_pos,i_pos-i_pipe,-1
	            element_counter=element_counter+1
                pipe_line(i,j_pos)=1.0
	            pipe_shape(i,j_pos:j_pos+j_radius)=1.0
	            pipe_shape(i,j_pos-j_radius:j_pos)=1.0
                l_pos=l_pos+dx_t
                l_line(element_counter)=l_pos
            end do
        end if
  
        if (pipe_counter>n_pipe) exit

        do j=j_pos,j_pos+j_pipe
            element_counter=element_counter+1
            pipe_line(i,j)=1.0
            pipe_shape(i:i+i_radius,j)=1.0
            pipe_shape(i-i_radius:i,j)=1.0
            l_pos=l_pos+dy_t
            l_line(element_counter)=l_pos
        end do
        i_pos=i
        j_pos=j_pos+j_pipe
    end do
    total_element=element_counter 
    i_marginf=i
    j_marginf=j

    allocate(kr(1:total_element),kl(1:total_element),al(1:total_element),ar(1:total_element),apt(1:total_element),ap0t(1:total_element),spt(1:total_element),sct(1:total_element),tt1(1:total_element),tf1(1:total_element),tf1_old(1:total_element),tt1_old(1:total_element),tf1_pre(1:total_element),tt1_pre(1:total_element),apf(1:total_element),ap0f(1:total_element),spf(1:total_element),scf(1:total_element),fr(1:total_element),fl(1:total_element),dr(1:total_element),dl(1:total_element),pr(1:total_element),pl(1:total_element),en1_t(1:total_element),en1_f(1:total_element))
    
    ! write shape of pipe to file 
    open(100,file='route.plt')
    write(100,*)'variables="x","y","Pipe Line","Pipe Shape"'
    write(100,*)'zone, i=', imax, 'j=', jmax
    do j=1,jmax
        do i=1,imax
            write(100,'(4f20.6)')x(i),y(j),pipe_line(i,j),pipe_shape(i,j)
        end do
    end do
    close(100)
    pause
    return
end subroutine
!-------------------------------------------------------------------------input properties---------------------------------------------------------------------------
subroutine input_property
    time=0.0
    sigma=5.670*(10**(-8.0))  !w/(m^2k^4)
    v_inf=1.0                   !m/s
    h_wind=h_wind_func(v_inf)
    tsky=0.0522*(tamb**1.5)   !k
    tsun=5800.0

    !************glass**************
    rho_g=2200.0         !kg/m^3
    c_g=480.0            !j/kg
    k_g=1.1           !w/mk
    alpha_g=0.05      
    eps_g=0.92
    tau_g=0.936  !1.0-alpha_g  
    d_g=0.003            !m thickness

    !*************air***************** 
    k_a=0.026               !w/mk
    v_a=16.0*(10**-6.0)     !m^2/s
    d_a=0.001                !m thickness
 
    !*************pv******************
    pac=0.94         !packing factor
    beta=0.0045       !% solar cell temp. coeff.
    eta0=16.0*10**-2.0        !% electrical conversion efficiency at reference temp. 298k
    rho_pv=2330.0    !kg/m^3
    c_pv=700.0       !j/kg
    k_pv=84.0       !w/mk
    alpha_pv=0.95      
    d_pv=0.0003       !m  thickness

    !*************adhesive************
    k_ad=200.0       !w/mk
    d_ad=0.0001       !m  thickness
    k1=g*tau_g*alpha_pv*eta0*pac
    k2=k1+k1*beta*298.0
    
    !********absorber plate*********** 
    rho_ab=8920.0     !kg/m^3
    c_ab=385.0
    d_ab=0.0004       !m thickness
    k_ab=398.0       !w/mk
    
    !***************tube**************
    rho_t=8920.0    !kg/m^3
    c_t=385.0       !j/kgk
    k_t=398.0       !w/mk
    per=pi*diameter      !m peripheral
    d_t=0.001      !m thickness
    a_t=per*d_t     !m^2 cross section area of tube
    a_t_ab=pi*diameter/8.0           ! contact area between tube and absorber per length
    a_t_in=pi*diameter*7.0/8.0       ! contact area between tube and insulation per length

    !*************fluid*************** 
    a_f=pi*(diameter**2.0)/4.0                             
    !v_f=0.1673

    !*************nano particles*************** 
    !al2o3
    rho_n=3970.0     !kg/m^3
    c_n=765.0         !j/kgk (nano)
    k_n=40.0           !w/mK
    d_n=2.0*10.0**(-8.0)

    !zno
    !rho_n=5600.0     !kg/m^3
    !c_n=495.0           !j/kgk (nano)
    !k_n=13.0           !w/mK
    !d_n=2.0*10.0**(-8.0)
    
    !tio2
    !rho_n=4250.0     !kg/m^3
    !c_n=686.0           !j/kgk (nano)
    !k_n=8.9           !w/mK
    !d_n=2.0*10.0**(-8.0)

    !Sio2
    !rho_n=2170.0     !kg/m^3
    !c_n=680.0           !j/kgk (nano)
    !k_n=1.3           !w/mK
    !d_n=2.0*10.0**(-8.0)

    !************insulation***********
    rho_i=48.0      !kg/m^3
    c_i=843.0       !j/kgk
    k_i=0.03749      !w/mk
    d_i=0.03        !m  thickness

    !*************conduction heat transfer***************************
    h_t_i=2.0/((d_i/k_i)+(d_t/k_t))
    h_ab_t=2.0/((d_ab/k_ab)+(d_t/k_t))
    h_ab_i=2.0/((d_ab/k_ab)+(d_i/k_i))
    h_pv_ab=2.0/((d_ab/k_ab)+(d_pv/k_pv)+(d_ad/k_ad))
    h_g_pv=0.3*2.0/((d_pv/k_pv)+(d_g/k_g)+(d_a/k_a))
    return
end subroutine
!----------------------------------------------------------------------------initialize------------------------------------------------------------------------------
subroutine initialize
    tg=tamb+0.01*g
    tpv=tamb+0.01*g
    tab=tamb+0.01*g
    ti=tamb+0.01*g
    tt1=tamb+0.01*g
    tt=tamb
    tf1=tf_inlet
    tf=tf_inlet

    tg_old=tg
    tpv_old=tpv
    tab_old=tab
    ti_old=ti
    tt1_old=tt1
    tf1_old=tf1

    tg_pre=tg
    tpv_pre=tpv
    tab_pre=tab
    ti_pre=ti
    tt1_pre=tt1
    tf1_pre=tf1
    return
end subroutine

subroutine fluid_property
    !*************fluid*************** 
    tave=sum(tf1)/size(tf1)
    mu_f=-0.000000003*tave**3.0+0.000003*tave**2.0-0.001*tave+0.1118
    k_f=-0.000009*(tave-273.0)**2.0+0.0021*(tave-273.0)+0.5603          !w/m^2k
    c_f=0.000005*(tave-273.0)**4.0-0.001*(tave-273.0)**3.0+0.0855*(tave-273.0)**2.0-3.0451*(tave-273.0)+4216
    rho_f=-0.0036*(tave-273.0)**2.0-0.0928*(tave-273.0)+1001.7     !kg/m^3

    !*************nanofluid*************** 
    phi=rho_f*phiw/rho_n
    mu_nf=mu_f*(1.0-phi)**(-2.5)
    rho_nf=phi*rho_n+(1.0-phi)*rho_f
    c_nf=(phi*rho_n*c_n+(1.0-phi)*rho_f*c_f)/(rho_nf)
    if(phi==0.0) then
        k_nf=k_f
    else
        !k_fx=(-0.8467*phi+0.0753)*tave+(237.67*phi-21.998)
        !k_beta=0.0017*(100*phi)**(-0.0841)
        !k_nf=k_f*(1.0+3.0*((k_n/k_f)-1.0)*phi/((k_n/k_f)+2.0-((k_n/k_f)-1.0)*phi))+50000*k_beta*phi*rho_n*c_n*k_fx*sqrt(1.38*10.0**(-23.0)*tave/(rho_n*d_n))
        k_nf=k_f*(k_n+2.0*k_f+2.0*(k_n-k_f)*phi)/(k_n+2.0*k_f-(k_n-k_f)*phi)
    end if
    v_f=4.0*m_f/(rho_nf*pi*diameter**2.0)
    !m_f=rho_nf*v_f*pi*diameter**2.0/4.0
    re_nf=4*m_f/(pi*diameter*mu_nf)
    pr_nf=mu_nf*c_nf/k_nf
    if(re_nf<3000) then
        nu_f=4.36+(0.086*(re_nf*pr_nf*diameter/l_line(total_element))**1.33)/(1.0+pr_nf*(re_nf*diameter/l_line(total_element))**0.83)
        !nu_f=4.36
    else
        f_nf=(0.79*log(re_nf)-1.64)**(-2.0)
        nu_f=((f_nf/8.0)*(re_nf-1000.0)*pr_nf)/(1.0+12.7*(f_nf/8.0)**0.5*(pr_nf**0.667-1.0))
    end if
    !*************convection*************** 
    hconvt_f=nu_f*k_nf/diameter
    
    return
end subroutine

subroutine open_files
    open(1,file=folder_array(file)//'_error_iteration.plt')
    open(2,file=folder_array(file)//'_error_time.plt')
    write(1,*)'variables="Iteration","Error"'
    write(2,*)'variables="Time Iteration","Error"'
    return
end subroutine

subroutine coefficients_glass
    do j=1,jmax
        do i=1,imax
            a=rho_g*c_g*dx(i)*dy(j)*d_g/dt
            b=dx(i)*dy(j)*(h_g_pv*tpv(i,j)+hrg_env(tg(i,j))*tsky+h_wind*tamb+alpha_g*g)
            c=dx(i)*dy(j)*(h_g_pv+hrg_env(tg(i,j))+h_wind)
            dn=k_g*dx(i)*d_g/dy_pn(j)
            ds=k_g*dx(i)*d_g/dy_sp(j)
            dw=k_g*dy(j)*d_g/dx_wp(i)
            de=k_g*dy(j)*d_g/dx_pe(i)
            if(i>1.and.i<imax.and.j>1.and.j<jmax) then !9
                aw(i,j)=-dw
                ae(i,j)=-de
                ap(i,j)=de+dw+dn+ds+a+c
	            sc(i,j)=b+a*tg_old(i,j)+dn*tg_old(i,j+1)+ds*tg_old(i,j-1)
            else if(i==1.and.j==1) then !3
                b=b+h_wind*tamb*(dx(i)+dy(j))*d_g
                c=c+h_wind*(dx(i)+dy(j))*d_g
                aw(i,j)=0
                ae(i,j)=-de
                ap(i,j)=de+dn+a+c
	            sc(i,j)=b+a*tg_old(i,j)+dn*tg_old(i,j+1)
            else if(i==1.and.j>1.and.j<jmax) then !2
                b=b+h_wind*tamb*dy(j)*d_g
                c=c+h_wind*dy(j)*d_g
                aw(i,j)=0
                ae(i,j)=-de
                ap(i,j)=de+dn+ds+a+c
	            sc(i,j)=b+a*tg_old(i,j)+dn*tg_old(i,j+1)+ds*tg_old(i,j-1)
            else if(i==1.and.j==jmax) then !1
                b=b+h_wind*tamb*(dx(i)+dy(j))*d_g
                c=c+h_wind*(dx(i)+dy(j))*d_g
                aw(i,j)=0
                ae(i,j)=-de
                ap(i,j)=de+ds+a+c
	            sc(i,j)=b+a*tg_old(i,j)+ds*tg_old(i,j-1)
            else if(i>1.and.i<imax.and.j==1) then !4
                b=b+h_wind*tamb*dx(i)*d_g
                c=c+h_wind*dx(i)*d_g
                aw(i,j)=-dw
                ae(i,j)=-de
                ap(i,j)=de+dw+dn+a+c
	            sc(i,j)=b+a*tg_old(i,j)+dn*tg_old(i,j+1)
            else if(i>1.and.i<imax.and.j==jmax) then !6
                b=b+h_wind*tamb*dy(j)*d_g
                c=c+h_wind*dy(j)*d_g
                aw(i,j)=-dw
                ae(i,j)=-de
                ap(i,j)=de+dw+ds+a+c
	            sc(i,j)=b+a*tg_old(i,j)+ds*tg_old(i,j-1)
            else if(i==imax.and.j==1) then !5
                b=b+h_wind*tamb*(dx(i)+dy(j))*d_g
                c=c+h_wind*(dx(i)+dy(j))*d_g
                aw(i,j)=-dw
                ae(i,j)=0
                ap(i,j)=dw+dn+a+c
	            sc(i,j)=b+a*tg_old(i,j)+dn*tg_old(i,j+1)
            else if(i==imax.and.j>1.and.j<jmax) then !8
                b=b+h_wind*tamb*dx(i)*d_g
                c=c+h_wind*dx(i)*d_g
                aw(i,j)=-dw
                ae(i,j)=0
                ap(i,j)=dw+dn+ds+a+c
	            sc(i,j)=b+a*tg_old(i,j)+dn*tg_old(i,j+1)+ds*tg_old(i,j-1)
            else if(i==imax.and.j==jmax) then !7
                b=b+h_wind*tamb*(dx(i)+dy(j))*d_g
                c=c+h_wind*(dx(i)+dy(j))*d_g
                aw(i,j)=-dw
                ae(i,j)=0
                ap(i,j)=dw+ds+a+c
	            sc(i,j)=b+a*tg_old(i,j)+ds*tg_old(i,j-1)
            end if
        end do
    end do
    return
end subroutine

subroutine coefficients_pv
    do j=1,jmax
        do i=1,imax
            a=rho_pv*c_pv*dx(i)*dy(j)*d_pv/dt
            b=dx(i)*dy(j)*(alpha_pv*g*tau_g+h_g_pv*tg(i,j)+h_pv_ab*tab(i,j)-k2)
            c=dx(i)*dy(j)*(h_g_pv+h_pv_ab-k1*beta)
            dw=k_pv*dy(j)*d_pv/dx_wp(i)
            de=k_pv*dy(j)*d_pv/dx_pe(i)
            dn=k_pv*dx(i)*d_pv/dy_pn(j)
            ds=k_pv*dx(i)*d_pv/dy_sp(j)
            if(i>1.and.i<imax.and.j>1.and.j<jmax) then !9
                aw(i,j)=-dw
                ae(i,j)=-de
                ap(i,j)=de+dw+dn+ds+a+c
	            sc(i,j)=b+a*tpv_old(i,j)+dn*tpv_old(i,j+1)+ds*tpv_old(i,j-1)
            else if(i==1.and.j==1) then !3
                b=b+h_wind*tamb*(dx(i)+dy(j))*d_pv
                c=c+h_wind*(dx(i)+dy(j))*d_pv
                aw(i,j)=0
                ae(i,j)=-de
                ap(i,j)=de+dn+a+c
	            sc(i,j)=b+a*tpv_old(i,j)+dn*tpv_old(i,j+1)
            else if(i==1.and.j>1.and.j<jmax) then !2
                b=b+h_wind*tamb*dy(j)*d_pv
                c=c+h_wind*dy(j)*d_pv
                aw(i,j)=0
                ae(i,j)=-de
                ap(i,j)=de+dn+ds+a+c
	            sc(i,j)=b+a*tpv_old(i,j)+dn*tpv_old(i,j+1)+ds*tpv_old(i,j-1)
            else if(i==1.and.j==jmax) then !1
                b=b+h_wind*tamb*(dx(i)+dy(j))*d_pv
                c=c+h_wind*(dx(i)+dy(j))*d_pv
                aw(i,j)=0
                ae(i,j)=-de
                ap(i,j)=de+ds+a+c
	            sc(i,j)=b+a*tpv_old(i,j)+ds*tpv_old(i,j-1)
            else if(i>1.and.i<imax.and.j==1) then !4
                b=b+h_wind*tamb*dx(i)*d_pv
                c=c+h_wind*dx(i)*d_pv
                aw(i,j)=-dw
                ae(i,j)=-de
                ap(i,j)=de+dw+dn+a+c
	            sc(i,j)=b+a*tpv_old(i,j)+dn*tpv_old(i,j+1)
            else if(i>1.and.i<imax.and.j==jmax) then !6
                b=b+h_wind*tamb*dy(j)*d_pv
                c=c+h_wind*dy(j)*d_pv
                aw(i,j)=-dw
                ae(i,j)=-de
                ap(i,j)=de+dw+ds+a+c
	            sc(i,j)=b+a*tpv_old(i,j)+ds*tpv_old(i,j-1)
            else if(i==imax.and.j==1) then !5
                b=b+h_wind*tamb*(dx(i)+dy(j))*d_pv
                c=c+h_wind*(dx(i)+dy(j))*d_pv
                aw(i,j)=-dw
                ae(i,j)=0
                ap(i,j)=dw+dn+a+c
	            sc(i,j)=b+a*tpv_old(i,j)+dn*tpv_old(i,j+1)
            else if(i==imax.and.j>1.and.j<jmax) then !8
                b=b+h_wind*tamb*dx(i)*d_pv
                c=c+h_wind*dx(i)*d_pv
                aw(i,j)=-dw
                ae(i,j)=0
                ap(i,j)=dw+dn+ds+a+c
	            sc(i,j)=b+a*tpv_old(i,j)+dn*tpv_old(i,j+1)+ds*tpv_old(i,j-1)
            else if(i==imax.and.j==jmax) then !7
                b=b+h_wind*tamb*(dx(i)+dy(j))*d_pv
                c=c+h_wind*(dx(i)+dy(j))*d_pv
                aw(i,j)=-dw
                ae(i,j)=0
                ap(i,j)=dw+ds+a+c
	            sc(i,j)=b+a*tpv_old(i,j)+ds*tpv_old(i,j-1)
            end if

        end do
    end do
    return
end subroutine

subroutine coefficients_absorber
    do j=1,jmax
        do i=1,imax
            a=rho_ab*c_ab*dx(i)*dy(j)*d_ab/dt
            b=dx(i)*dy(j)*(h_pv_ab*tpv(i,j)+pipe_line(i,j)*(a_t_ab/dy_t)*h_ab_t*tt(i,j)+(1.0-pipe_line(i,j))*h_ab_i*ti(i,j))
            c=dx(i)*dy(j)*(h_pv_ab+pipe_line(i,j)*(a_t_ab/dy_t)*h_ab_t+(1.0-pipe_line(i,j))*h_ab_i)
            dw=k_ab*dy(j)*d_ab/dx_wp(i)
            de=k_ab*dy(j)*d_ab/dx_pe(i)
            dn=k_ab*dx(i)*d_ab/dy_pn(j)
            ds=k_ab*dx(i)*d_ab/dy_sp(j)
            if(i>1.and.i<imax.and.j>1.and.j<jmax) then !9
                aw(i,j)=-dw
                ae(i,j)=-de
                ap(i,j)=de+dw+dn+ds+a+c
	            sc(i,j)=b+a*tab_old(i,j)+dn*tab_old(i,j+1)+ds*tab_old(i,j-1)
            else if(i==1.and.j==1) then !3
                b=b+h_wind*tamb*(dx(i)+dy(j))*d_ab
                c=c+h_wind*(dx(i)+dy(j))*d_ab
                aw(i,j)=0
                ae(i,j)=-de
                ap(i,j)=de+dn+a+c
	            sc(i,j)=b+a*tab_old(i,j)+dn*tab_old(i,j+1)
            else if(i==1.and.j>1.and.j<jmax) then !2
                b=b+h_wind*tamb*dy(j)*d_ab
                c=c+h_wind*dy(j)*d_ab
                aw(i,j)=0
                ae(i,j)=-de
                ap(i,j)=de+dn+ds+a+c
	            sc(i,j)=b+a*tab_old(i,j)+dn*tab_old(i,j+1)+ds*tab_old(i,j-1)
            else if(i==1.and.j==jmax) then !1
                b=b+h_wind*tamb*(dx(i)+dy(j))*d_ab
                c=c+h_wind*(dx(i)+dy(j))*d_ab
                aw(i,j)=0
                ae(i,j)=-de
                ap(i,j)=de+ds+a+c
	            sc(i,j)=b+a*tab_old(i,j)+ds*tab_old(i,j-1)
            else if(i>1.and.i<imax.and.j==1) then !4
                b=b+h_wind*tamb*dx(i)*d_ab
                c=c+h_wind*dx(i)*d_ab
                aw(i,j)=-dw
                ae(i,j)=-de
                ap(i,j)=de+dw+dn+a+c
	            sc(i,j)=b+a*tab_old(i,j)+dn*tab_old(i,j+1)
            else if(i>1.and.i<imax.and.j==jmax) then !6
                b=b+h_wind*tamb*dy(j)*d_ab
                c=c+h_wind*dy(j)*d_ab
                aw(i,j)=-dw
                ae(i,j)=-de
                ap(i,j)=de+dw+ds+a+c
	            sc(i,j)=b+a*tab_old(i,j)+ds*tab_old(i,j-1)
            else if(i==imax.and.j==1) then !5
                b=b+h_wind*tamb*(dx(i)+dy(j))*d_ab
                c=c+h_wind*(dx(i)+dy(j))*d_ab
                aw(i,j)=-dw
                ae(i,j)=0
                ap(i,j)=dw+dn+a+c
	            sc(i,j)=b+a*tab_old(i,j)+dn*tab_old(i,j+1)
            else if(i==imax.and.j>1.and.j<jmax) then !8
                b=b+h_wind*tamb*dx(i)*d_ab
                c=c+h_wind*dx(i)*d_ab
                aw(i,j)=-dw
                ae(i,j)=0
                ap(i,j)=dw+dn+ds+a+c
	            sc(i,j)=b+a*tab_old(i,j)+dn*tab_old(i,j+1)+ds*tab_old(i,j-1)
            else if(i==imax.and.j==jmax) then !7
                b=b+h_wind*tamb*(dx(i)+dy(j))*d_ab
                c=c+h_wind*(dx(i)+dy(j))*d_ab
                aw(i,j)=-dw
                ae(i,j)=0
                ap(i,j)=dw+ds+a+c
	            sc(i,j)=b+a*tab_old(i,j)+ds*tab_old(i,j-1)
            end if
        end do
    end do
    return
end subroutine

subroutine coefficients_tube
    i_pos=i_margin
    j_pos=j_margin
    pipe_counter=1
    element_counter=0
    do
        pipe_counter=pipe_counter+1
        if(pipe_counter==9) then  
            do i=i_pos,i_pos-i_center,-1
                element_counter=element_counter+1
	            dv_t=a_t*dx(i)
                a=rho_t*c_t*dv_t/dt
                b=(a_t_ab*h_ab_t*tab(i,j_pos)+per*hconvt_f*tf(i,j_pos)+a_t_in*h_t_i*ti(i,j_pos))*dx(i)
                c=(a_t_ab*h_ab_t+per*hconvt_f+a_t_in*h_t_i)*dx(i)
                if(element_counter==total_element) then
                    de=0
                else
                    de=k_t*a_t/dx_pe(i)
                end if
                dw=k_t*a_t/dx_wp(i)
	            ar(element_counter)=-de
	            al(element_counter)=-dw
                apt(element_counter)=de+dw+c+a
                sct(element_counter)=b+a*tt1_old(element_counter)
                
            end do
        else if(pipe_counter==10) then  
            do i=i_pos,i_pos+i_center,1
                element_counter=element_counter+1
                dv_t=a_t*dx(i)
                a=rho_t*c_t*dv_t/dt
                b=(a_t_ab*h_ab_t*tab(i,j_pos)+per*hconvt_f*tf(i,j_pos)+a_t_in*h_t_i*ti(i,j_pos))*dx(i)
                c=(a_t_ab*h_ab_t+per*hconvt_f+a_t_in*h_t_i)*dx(i)
                de=k_t*a_t/dx_pe(i)
                if(element_counter==1) then
                    dw=0
                else
                    dw=k_t*a_t/dx_wp(i)
                end if
	            ar(element_counter)=-de
	            al(element_counter)=-dw
                apt(element_counter)=de+dw+c+a
                sct(element_counter)=b+a*tt1_old(element_counter)
            end do
        else if(mod(pipe_counter,2)==0) then  
            do i=i_pos,i_pos+i_pipe,1
                element_counter=element_counter+1
                dv_t=a_t*dx(i)
                a=rho_t*c_t*dv_t/dt
                b=(a_t_ab*h_ab_t*tab(i,j_pos)+per*hconvt_f*tf(i,j_pos)+a_t_in*h_t_i*ti(i,j_pos))*dx(i)
                c=(a_t_ab*h_ab_t+per*hconvt_f+a_t_in*h_t_i)*dx(i)
                de=k_t*a_t/dx_pe(i)
                if(element_counter==1) then
                    dw=0
                else
                    dw=k_t*a_t/dx_wp(i)
                end if
	            ar(element_counter)=-de
	            al(element_counter)=-dw
                apt(element_counter)=de+dw+c+a
                sct(element_counter)=b+a*tt1_old(element_counter)
            end do
        else
            do i=i_pos,i_pos-i_pipe,-1
                element_counter=element_counter+1
	            dv_t=a_t*dx(i)
                a=rho_t*c_t*dv_t/dt
                b=(a_t_ab*h_ab_t*tab(i,j_pos)+per*hconvt_f*tf(i,j_pos)+a_t_in*h_t_i*ti(i,j_pos))*dx(i)
                c=(a_t_ab*h_ab_t+per*hconvt_f+a_t_in*h_t_i)*dx(i)
                if(element_counter==total_element) then
                    de=0
                else
                    de=k_t*a_t/dx_pe(i)
                end if
                dw=k_t*a_t/dx_wp(i)
	            ar(element_counter)=-de
	            al(element_counter)=-dw
                apt(element_counter)=de+dw+c+a
                sct(element_counter)=b+a*tt1_old(element_counter)
                
            end do
        end if
  
        if (pipe_counter>n_pipe) exit

        do j=j_pos,j_pos+j_pipe
            element_counter=element_counter+1
            dv_t=a_t*dy(j)
            a=rho_t*c_t*dv_t/dt
            b=(a_t_ab*h_ab_t*tab(i,j)+per*hconvt_f*tf(i,j)+a_t_in*h_t_i*ti(i,j))*dy(j)
            c=(a_t_ab*h_ab_t+per*hconvt_f+a_t_in*h_t_i)*dy(j)
            de=k_t*a_t/dy_pn(j)
            dw=k_t*a_t/dy_sp(j)
            ar(element_counter)=-de
            al(element_counter)=-dw
            apt(element_counter)=de+dw+c+a
            sct(element_counter)=b+a*tt1_old(element_counter)
        end do
        i_pos=i
        j_pos=j_pos+j_pipe
    end do  	 
    return
end subroutine

subroutine res_tube

    !--------------------assign the tt1(:) to tt(:,:) by considering the rout------------------
    i_pos=i_margin
    j_pos=j_margin
    pipe_counter=1
    element_counter=0
    do
        pipe_counter=pipe_counter+1
        if(pipe_counter==9) then  
            do i=i_pos,i_pos-i_center,-1
                element_counter=element_counter+1
                tt(i,j_pos)=tt1(element_counter)	  
            end do
        else if(pipe_counter==10) then  
            do i=i_pos,i_pos+i_center,1
	            element_counter=element_counter+1
                tt(i,j_pos)=tt1(element_counter)
	        end do
        else if(mod(pipe_counter,2)==0) then  
            do i=i_pos,i_pos+i_pipe,1
	            element_counter=element_counter+1
                tt(i,j_pos)=tt1(element_counter)
	        end do
        else
            do i=i_pos,i_pos-i_pipe,-1
                element_counter=element_counter+1
                tt(i,j_pos)=tt1(element_counter)	  
            end do
        end if
  
        if (pipe_counter>n_pipe) exit

        do j=j_pos,j_pos+j_pipe
            element_counter=element_counter+1
            tt(i,j)=tt1(element_counter)    
        end do
        i_pos=i
        j_pos=j_pos+j_pipe
    end do  	 
    return
end subroutine

subroutine coefficients_fluid
    gamma=k_nf/c_nf
    i_pos=i_margin
    j_pos=j_margin
    pipe_counter=1
    element_counter=0
    do
        pipe_counter=pipe_counter+1
        if(pipe_counter==9) then
            do i=i_pos,i_pos-i_center,-1
                element_counter=element_counter+1
	            dv_f=a_f*dx(i)
                a=rho_nf*c_nf*dv_f/dt
                b=1.0*per*hconvt_f*dx(i)*tt1(element_counter)
                c=1.0*per*hconvt_f*dx(i)+m_f*c_nf
                if(element_counter==total_element) then
                    de=0
                else
                    de=k_nf*a_f/dx_pe(i)
                end if
                dw=k_nf*a_f/dx_wp(i)
                ar(element_counter)=-de
                al(element_counter)=-dw-m_f*c_nf
                apf(element_counter)=a+c+de+dw
                scf(element_counter)=b+a*tf1_old(element_counter)
            end do
        else if(pipe_counter==10) then
            do i=i_pos,i_pos+i_center,1
                element_counter=element_counter+1
	            dv_f=a_f*dx(i)
                a=rho_nf*c_nf*dv_f/dt
                b=1.0*per*hconvt_f*dx(i)*tt1(element_counter)
                c=1.0*per*hconvt_f*dx(i)+m_f*c_nf
                de=k_nf*a_f/dx_pe(i)
                dw=k_nf*a_f/dx_wp(i)
                if(element_counter==1) then
                    b=b+m_f*c_nf*tf_inlet+dw*tf_inlet
                    ar(element_counter)=-de
                    al(element_counter)=0
                    apf(element_counter)=a+c+de+dw
                    scf(element_counter)=b+a*tf1_old(element_counter)
                else
                    ar(element_counter)=-de
                    al(element_counter)=-dw-m_f*c_nf
                    apf(element_counter)=a+c+de+dw
                    scf(element_counter)=b+a*tf1_old(element_counter)                    
                end if
            end do
        else if(mod(pipe_counter,2)==0) then
            do i=i_pos,i_pos+i_pipe,1
                element_counter=element_counter+1
	            dv_f=a_f*dx(i)
                a=rho_nf*c_nf*dv_f/dt
                b=1.0*per*hconvt_f*dx(i)*tt1(element_counter)
                c=1.0*per*hconvt_f*dx(i)+m_f*c_nf
                de=k_nf*a_f/dx_pe(i)
                dw=k_nf*a_f/dx_wp(i)
                if(element_counter==1) then
                    b=b+m_f*c_nf*tf_inlet+dw*tf_inlet
                    ar(element_counter)=-de
                    al(element_counter)=0
                    apf(element_counter)=a+c+de+dw
                    scf(element_counter)=b+a*tf1_old(element_counter)
                else
                    ar(element_counter)=-de
                    al(element_counter)=-dw-m_f*c_nf
                    apf(element_counter)=a+c+de+dw
                    scf(element_counter)=b+a*tf1_old(element_counter)                    
                end if
            end do
        else
            do i=i_pos,i_pos-i_pipe,-1
                element_counter=element_counter+1
	            dv_f=a_f*dx(i)
                a=rho_nf*c_nf*dv_f/dt
                b=1.0*per*hconvt_f*dx(i)*tt1(element_counter)
                c=1.0*per*hconvt_f*dx(i)+m_f*c_nf
                if(element_counter==total_element) then
                    de=0
                else
                    de=k_nf*a_f/dx_pe(i)
                end if
                dw=k_nf*a_f/dx_wp(i)
                ar(element_counter)=-de
                al(element_counter)=-dw-m_f*c_nf
                apf(element_counter)=a+c+de+dw
                scf(element_counter)=b+a*tf1_old(element_counter)
            end do
        end if
  
        if (pipe_counter>n_pipe) exit

        do j=j_pos,j_pos+j_pipe
            element_counter=element_counter+1
            dv_f=a_f*dy(j)
            a=rho_nf*c_nf*dv_f/dt
            b=1.0*per*hconvt_f*dy(j)*tt1(element_counter)
            c=1.0*per*hconvt_f*dy(j)+m_f*c_nf
            de=k_nf*a_f/dy_pn(j)
            dw=k_nf*a_f/dy_sp(j)
            ar(element_counter)=-de
            al(element_counter)=-dw-m_f*c_nf
            apf(element_counter)=a+c+de+dw
            scf(element_counter)=b+a*tf1_old(element_counter)
        end do
        i_pos=i
        j_pos=j_pos+j_pipe
    end do  	 
    return
end subroutine

subroutine res_fluid

    !--------------------assign the tf1(:) to tf(:,:) by considering the rout------------------
    i_pos=i_margin
    j_pos=j_margin
    pipe_counter=1
    element_counter=0
    do
        pipe_counter=pipe_counter+1
        if(pipe_counter==9) then  
            do i=i_pos,i_pos-i_center,-1
                element_counter=element_counter+1
                tf(i,j_pos)=tf1(element_counter)	  
            end do
        else if(pipe_counter==10) then  
            do i=i_pos,i_pos+i_center,1
	            element_counter=element_counter+1
                tf(i,j_pos)=tf1(element_counter)
	        end do
        else if(mod(pipe_counter,2)==0) then  
            do i=i_pos,i_pos+i_pipe,1
	            element_counter=element_counter+1
                tf(i,j_pos)=tf1(element_counter)
	        end do
        else
            do i=i_pos,i_pos-i_pipe,-1
                element_counter=element_counter+1
                tf(i,j_pos)=tf1(element_counter)	  
            end do
        end if
  
        if (pipe_counter>n_pipe) exit

        do j=j_pos,j_pos+j_pipe
            element_counter=element_counter+1
            tf(i,j)=tf1(element_counter)
        
        end do
        i_pos=i
        j_pos=j_pos+j_pipe
    end do  	 
return
end subroutine

subroutine coefficients_insulation
    do j=1,jmax
        do i=1,imax
            a=rho_i*c_i*dv(i,j)*d_i/dt
            b=(pipe_line(i,j)*(a_t_in/dy_t)*h_t_i*tt(i,j)+(1.0-pipe_line(i,j))*h_ab_i*tab(i,j)+h_wind*tamb)*dv(i,j)
            c=(pipe_line(i,j)*(a_t_in/dy_t)*h_t_i+(1.0-pipe_line(i,j))*h_ab_i+h_wind)*dv(i,j)
            dw=k_i*dy(j)*d_i/dx_wp(i)
            de=k_i*dy(j)*d_i/dx_pe(i)
            dn=k_i*dx(i)*d_i/dy_pn(j)
            ds=k_i*dx(i)*d_i/dy_sp(j)
            if(i>1.and.i<imax.and.j>1.and.j<jmax) then !9
                aw(i,j)=-dw
                ae(i,j)=-de
                ap(i,j)=de+dw+dn+ds+a+c
	            sc(i,j)=b+a*ti_old(i,j)+dn*ti_old(i,j+1)+ds*ti_old(i,j-1)
            else if(i==1.and.j==1) then !3
                b=b+h_wind*tamb*(dx(i)+dy(j))*d_i
                c=c+h_wind*(dx(i)+dy(j))*d_i
                aw(i,j)=0
                ae(i,j)=-de
                ap(i,j)=de+dn+a+c
	            sc(i,j)=b+a*ti_old(i,j)+dn*ti_old(i,j+1)
            else if(i==1.and.j>1.and.j<jmax) then !2
                b=b+h_wind*tamb*dy(j)*d_i
                c=c+h_wind*dy(j)*d_i
                aw(i,j)=0
                ae(i,j)=-de
                ap(i,j)=de+dn+ds+a+c
	            sc(i,j)=b+a*ti_old(i,j)+dn*ti_old(i,j+1)+ds*ti_old(i,j-1)
            else if(i==1.and.j==jmax) then !1
                b=b+h_wind*tamb*(dx(i)+dy(j))*d_i
                c=c+h_wind*(dx(i)+dy(j))*d_i
                aw(i,j)=0
                ae(i,j)=-de
                ap(i,j)=de+ds+a+c
	            sc(i,j)=b+a*ti_old(i,j)+ds*ti_old(i,j-1)
            else if(i>1.and.i<imax.and.j==1) then !4
                b=b+h_wind*tamb*dx(i)*d_i
                c=c+h_wind*dx(i)*d_i
                aw(i,j)=-dw
                ae(i,j)=-de
                ap(i,j)=de+dw+dn+a+c
	            sc(i,j)=b+a*ti_old(i,j)+dn*ti_old(i,j+1)
            else if(i>1.and.i<imax.and.j==jmax) then !6
                b=b+h_wind*tamb*dy(j)*d_i
                c=c+h_wind*dy(j)*d_i
                aw(i,j)=-dw
                ae(i,j)=-de
                ap(i,j)=de+dw+ds+a+c
	            sc(i,j)=b+a*ti_old(i,j)+ds*ti_old(i,j-1)
            else if(i==imax.and.j==1) then !5
                b=b+h_wind*tamb*(dx(i)+dy(j))*d_i
                c=c+h_wind*(dx(i)+dy(j))*d_i
                aw(i,j)=-dw
                ae(i,j)=0
                ap(i,j)=dw+dn+a+c
	            sc(i,j)=b+a*ti_old(i,j)+dn*ti_old(i,j+1)
            else if(i==imax.and.j>1.and.j<jmax) then !8
                b=b+h_wind*tamb*dx(i)*d_i
                c=c+h_wind*dx(i)*d_i
                aw(i,j)=-dw
                ae(i,j)=0
                ap(i,j)=dw+dn+ds+a+c
	            sc(i,j)=b+a*ti_old(i,j)+dn*ti_old(i,j+1)+ds*ti_old(i,j-1)
            else if(i==imax.and.j==jmax) then !7
                b=b+h_wind*tamb*(dx(i)+dy(j))*d_i
                c=c+h_wind*(dx(i)+dy(j))*d_i
                aw(i,j)=-dw
                ae(i,j)=0
                ap(i,j)=dw+ds+a+c
	            sc(i,j)=b+a*ti_old(i,j)+ds*ti_old(i,j-1)
            end if
        end do
    end do
    return
end subroutine

subroutine efficiency
    e_th=m_f*c_nf*(tf1(total_element)-tf1(1))/(l*h)
    e_in=tau_g*alpha_pv*g
    e_el_pv=-k1*beta*sum(tpv)/size(tpv)+k2
    e_el=e_el_pv
    etha_th=e_th/g*100.0
    etha_el_pv=e_el_pv/g*100.0
    etha_el=e_el/g*100.0
    etha_pvt=etha_th+etha_el
    etha_tht=etha_th+etha_el/0.38
    ex_in=g*(1.0-4.0*tamb/(3.0*tsun)+(tamb/tsun)**4.0/3.0)
    !ex_in=g*(1.0-tamb/tsun)
    ex_th=e_th-m_f*c_nf*tamb*log(tf1(total_element)/tf1(1))/(l*h)
    ex_el=e_el
    eps_pvt=(ex_th+ex_el)/ex_in*100.0
    eps_th=ex_th/ex_in*100.0
    eps_el=ex_el/ex_in*100.0
    en_system=(ex_in-ex_th-ex_el)/tamb
    return
end subroutine

subroutine entropy_fluid
    en_f_t=0.0
    en_f_f=0.0
    en1_t=0.0
    en1_f=0.0
    en1_t(1)=k_nf*(((tf1(1)-tf_inlet)/(l_line(1)))**2.0+(((tt1(1)+tf1(1))/2.0-tf1(1))/(diameter/2.0))**2.0)/tamb**2.0
    dp=128.0*m_f*mu_nf*(l_line(total_element)-l_line(1))/(rho_nf*pi*diameter**4.0)
    do k=2,total_element
        en1_t(k)=k_nf*(((tf1(k)-tf1(k-1))/(l_line(k)-l_line(k-1)))**2.0+(((tt1(k)+tf1(k))/2.0-tf1(k))/(diameter/2.0))**2.0)/tamb**2.0
    end do
    i_pos=i_margin
    j_pos=j_margin
    pipe_counter=1
    element_counter=0
    do
        pipe_counter=pipe_counter+1
        if(pipe_counter==9) then  
            do i=i_pos,i_pos-i_center,-1
                element_counter=element_counter+1
                en_f_t(i,j_pos)=en1_t(element_counter)
            end do
        else if(pipe_counter==10) then  
            do i=i_pos,i_pos+i_center,1
	            element_counter=element_counter+1
                en_f_t(i,j_pos)=en1_t(element_counter)	  
	        end do
        else if(mod(pipe_counter,2)==0) then  
            do i=i_pos,i_pos+i_pipe,1
	            element_counter=element_counter+1
                en_f_t(i,j_pos)=en1_t(element_counter)	  
	        end do
        else
            do i=i_pos,i_pos-i_pipe,-1
                element_counter=element_counter+1
                en_f_t(i,j_pos)=en1_t(element_counter)	  
            end do
        end if
  
        if (pipe_counter>n_pipe) exit

        do j=j_pos,j_pos+j_pipe
            element_counter=element_counter+1
            en_f_t(i,j)=en1_t(element_counter)
        end do
        i_pos=i
        j_pos=j_pos+j_pipe
    end do
    en_f_t_total=en1_t(1)*(l_line(1))*pi*diameter**2.0/4.0
    do k=2,total_element
        en_f_t_total=en_f_t_total+en1_t(k)*(l_line(k)-l_line(k-1))*pi*diameter**2.0/4.0
    end do
    en_f_t2_total=(m_f*c_nf*(tf1(total_element)-tf_inlet))**2.0/(k_nf*tamb**2.0*nu_f*l_line(total_element)**2.0)
    !en_f_f_total=8.0*m_f**3.0*64.0/(pi**2.0*rho_nf**2.0*tamb*diameter**5.0*re_nf)
    en_f_f_total=m_f*dp*log(tf1(total_element)/tf_inlet)/(rho_nf*(tf1(total_element)-tf_inlet))
    return
end subroutine

subroutine energy_balance
    p_s = g*l*h
    p_s_g = alpha_g*g*l*h
    p_s_pv = g*tau_g*l*h
    p_g_ref = (1.0-alpha_g-tau_g)*g*l*h
    p_pv = e_el_pv*l*h
    p_fluid = e_th*l*h
    p_g_conv = 0.0
    p_g_pv = 0.0
    p_g_sky = 0.0
    p_pv_ab = 0.0
    p_ab_t = 0.0
    p_ab_in = 0.0
    p_t_in = 0.0
    p_in_conv = 0.0
    do j=1,jmax
        do i=1,imax
            p_g_conv = p_g_conv + h_wind*(tg(i,j)-tamb)*dv(i,j)
            p_g_pv = p_g_pv + h_g_pv*(tg(i,j)-tpv(i,j))*dv(i,j)
            p_g_sky = p_g_sky + hrg_env(tg(i,j))*(tg(i,j)-tsky)*dv(i,j)
            p_pv_ab = p_pv_ab + h_pv_ab*(tpv(i,j)-tab(i,j))*dv(i,j)
            p_ab_t = p_ab_t + pipe_line(i,j)*(a_t_ab/dy_t)*h_ab_t*(tab(i,j)-tt(i,j))*dv(i,j)
            p_ab_in = p_ab_in + (1.0-pipe_line(i,j))*h_ab_i*(tab(i,j)-ti(i,j))*dv(i,j)
            p_t_in = p_t_in + pipe_line(i,j)*(a_t_in/dy_t)*h_t_i*(tt(i,j)-ti(i,j))*dv(i,j)
            p_in_conv = p_in_conv + h_wind*(ti(i,j)-tamb)*dv(i,j)
        end do
    end do
    
    return
end subroutine

subroutine output
    call cpu_time(end_time)
    call date_and_time(DATE=date)
    open(101,file=folder_array(file)//'_contour_glass.plt')
    open(102,file=folder_array(file)//'_contour_pv.plt')
    open(103,file=folder_array(file)//'_contour_absorber.plt')
    open(104,file=folder_array(file)//'_contour_tube.plt')
    open(105,file=folder_array(file)//'_contour_fluid.plt')
    open(106,file=folder_array(file)//'_contour_insulation.plt')
    open(107,file=folder_array(file)//'_contour_entropy_thermal_fluid.plt')
    write(101,*)'variables="x","y","Glass Temperature (K)"'
    write(102,*)'variables="x","y","PV Temperature (K)"'
    write(103,*)'variables="x","y","Absorber Temperature (K)"'
    write(104,*)'variables="x","y","Tube Temperature (K)"'
    write(105,*)'variables="x","y","Fluid Temperature (K)"'
    write(106,*)'variables="x","y","Insulation Temperature (K)"'
    write(107,*)'variables="x","y","Fluid Thermal Entropy Generation (w/Km^3)"'
    write(101,*)'zone, i=', imax, 'j=', jmax
    write(102,*)'zone, i=', imax, 'j=', jmax
    write(103,*)'zone, i=', imax, 'j=', jmax
    write(104,*)'zone, i=', imax, 'j=', jmax
    write(105,*)'zone, i=', imax, 'j=', jmax
    write(106,*)'zone, i=', imax, 'j=', jmax
    write(107,*)'zone, i=', imax, 'j=', jmax
    do j=1,jmax
        do i=1,imax
            write(101,'(3f50.20)')x(i),y(j),tg(i,j)
            write(102,'(3f50.20)')x(i),y(j),tpv(i,j)
            write(103,'(3f50.20)')x(i),y(j),tab(i,j)
            write(104,'(3f50.20)')x(i),y(j),tt(i,j)
            write(105,'(3f50.20)')x(i),y(j),tf(i,j)
            write(106,'(3f50.20)')x(i),y(j),ti(i,j)
            write(107,'(3f50.20)')x(i),y(j),en_f_t(i,j)
        end do
    end do
    write(200,'(a11,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,a1,f50.20,a1,f50.20,a1,f50.20,a1,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,a1,f50.20,a1,f50.20,a1,f50.20,a1,a1,a8,a1,f50.20,a1,f50.20,a1,f50.20,a1,i20,a1,i20)') folder_array(file),',',test_time(file),',',g,',',tamb,',',tf_inlet,',',phiw*100.0,',',phi*100.0,',',',',tf1(total_element),',',sum(tg)/size(tg),',',sum(tpv)/size(tpv),',',sum(tte)/size(tte),',',sum(tab)/size(tab),',',sum(tt1)/size(tt1),',',sum(tf1)/size(tf1),',',sum(ti)/size(ti),',',sum(tpv)/size(tpv)-sum(tab)/size(tab),',',',',maxval(tg),',',maxval(tpv),',',maxval(tte),',',maxval(tab),',',maxval(tt1),',',maxval(tf1),',',maxval(ti),',',',',ex_outlet(file),',',ex_surface(file),',',tf1(total_element)-ex_outlet(file),',',ex_surface(file)-maxval(tpv),',',',',etha_pvt,',',etha_th,',',etha_el_pv,',',etha_el_te,',',etha_el,',',ex_elef(file)*100.0,',',ex_elef(file)*100.0-etha_el,',',etha_tht,',',',',eps_th,',',eps_el,',',ex_th,',',ex_in,',',ex_el,',',en_f_t_total,',',en_f_t2_total,',',dp,',',en_f_f_total,',',en_system,',',',',re_nf,',',pr_nf,',',nu_f,',',',',hconvt_f,',',h_t_i,',',h_ab_t,',',h_ab_i,',',h_pv_te,',',h_te_ab,',',h_g_pv,',',',',g*l*h,',',alpha_g*g*l*h,',',h_g_pv*(sum(tg)/size(tg)-sum(tpv)/size(tpv))*l*h,',',g*tau_g*l*h,',',e_el_pv*l*h,',',h_pv_te*(sum(tpv)/size(tpv)-sum(tte)/size(tte))*l*h,',',e_el_te*l*h,',',h_te_ab*(sum(tte)/size(tte)-sum(tab)/size(tab))*l*h,',',h_ab_t*(sum(tab)/size(tab)-sum(tt1)/size(tt1))*l*h,',',e_th*l*h,',',',',rho_nf,',',c_nf,',',k_nf,',',mu_nf,',',',',tsky,',',v_f,',',m_f,',',',',date,',',end_time-begin_time,',',dt,',',time,',',imax,',',jmax
    write(201,'(26f50.20)') phi*100.0,phiw*100.0,test_time(file),tf1(total_element),maxval(tpv),ex_outlet(file),ex_surface(file),etha_pvt,etha_th,etha_el,ex_elef(file)*100.0,etha_tht,eps_pvt,eps_th,eps_el,rho_nf,c_nf,k_nf,mu_nf,re_nf,pr_nf,nu_f,en_f_t_total,en_f_t2_total,en_f_f_total,en_system
    write(202,'(f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20)') g,',',tf_inlet,',',m_f,',',tf1(total_element),',',sum(tpv)/size(tpv)-sum(tab)/size(tab),',',etha_pvt,',',etha_th,',',etha_el_pv,',',etha_el_te,',',etha_el,',',eps_th,',',e_el_pv*l*h,',',h_pv_te*(sum(tpv)/size(tpv)-sum(tte)/size(tte))*l*h,',',e_el_te*l*h,',',e_th*l*h
    write(203,'(f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20,a1,f50.20)') p_s,',',p_s_g,',',p_s_pv,',',p_g_ref,',',p_pv,',',p_fluid,',',p_g_conv,',',p_g_pv,',',p_g_sky,',',p_pv_ab,',',p_ab_t,',',p_ab_in,',',p_t_in,',',p_in_conv

    !close(101)
    !close(102)
    !close(103)
    !close(104)
    !close(105)
    !close(106)
    !close(107)
    close(1)
    close(2)
    return
end subroutine

subroutine tdma(a,b,c,d,n,x)
    implicit none
    double precision :: q
    double precision, dimension(n) :: x,a,b,c,d
    integer :: i,n
    do i=2,n
        q=a(i)/b(i-1)
        b(i)=b(i)-c(i-1)*q
        d(i)=d(i)-d(i-1)*q
    end do
    q=d(n)/b(n)
    x(n)=q
    do i=n-1,1,-1
        q=(d(i)-c(i)*q)/b(i)
        x(i)=q
    end do
end subroutine

double precision function h_wind_func(xx)
    double precision,intent(in) ::xx
    if (xx<=5.0) then
        h_wind_func=5.70+3.8*xx
    else
        h_wind_func=6.47+xx**0.78
    end if
end function

double precision function hrg_env(xx)
    double precision,intent(in) ::xx
    hrg_env=eps_g*sigma*(xx**2.0+tsky**2.0)*(xx+tsky)
end function

double precision function ah_d(xx)
    double precision,intent(in) ::xx
    ah_d=max(0.0,1.0-0.5*abs(xx))
end function

double precision function x_func(i)
    integer,intent(in) :: i
    if(i==1) then
        x_func=0.0
    else if(i==2) then
        x_func=dx_t/2.0
    else if(i>2 .and. i<=(imax-1)) then
        x_func=dx_t/2.0+(i-2.0)*dx_t
    else
        x_func=(i-1.0)*dx_t
    end if
end function

double precision function y_func(j)
    integer,intent(in) :: j
    if(j==1) then
        y_func=0.0
    else if(j==2) then
        y_func=dy_t/2.0
    else if(j>2 .and. j<=(jmax-1)) then
        y_func=dy_t/2.0+(j-2.0)*dy_t
    else
        y_func=(j-1.0)*dy_t
    end if
end function

double precision function dx_pe_func(i)
    integer,intent(in) :: i
    dx_pe_func=x(i+1)-x(i)
end function

double precision function dx_wp_func(i)
    integer,intent(in) :: i
    dx_wp_func=x(i)-x(i-1)
end function

double precision function dy_pn_func(j)
    integer,intent(in) :: j
    dy_pn_func=y(j+1)-y(j)
end function

double precision function dy_sp_func(j)
    integer,intent(in) :: j
    dy_sp_func=y(j)-y(j-1)
end function

double precision function dy_func(j)
    integer,intent(in) :: j
    if(j==1) then
        dy_func=dy_t/4.0
    else if(j==2) then
        dy_func=dy_sp(2)+dy_pn(2)/2.0
    else if(j>2 .and. j<(jmax-1)) then
        dy_func=(y(j+1)-y(j-1))/2.0
    else if(j==jmax-1) then
        dy_func=dy_pn(jmax-1)+dy_sp(jmax-1)/2.0
    else
        dy_func=dy_t/4.0
    end if
end function

double precision function dx_func(i)
    integer,intent(in) :: i
    if(i==1) then
        dx_func=dx_t/4.0
    else if(i==2) then
        dx_func=dx_wp(2)+dx_pe(2)/2.0
    else if(i>2 .and. i<(imax-1)) then
        dx_func=(x(i+1)-x(i-1))/2.0
    else if(i==imax-1) then
        dx_func=dx_pe(imax-1)+dx_wp(imax-1)/2.0
    else
        dx_func=dx_t/4.0
    end if
end function

double precision function dv_func(i,j)
    integer,intent(in) :: i,j
    dv_func=dy(j)*dx(i)
end function

end program pvt