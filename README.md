# mosel_test

 model ModelName
uses "mmxprs", "mmsystem"; !gain access to the Xpress-Optimizer solver

	!Parameters
	parameters
		DIR='./data'
		FILETYPE='.dat'
		datasizemax=20
		DELIM=','
    	D_TYPE='separate'
    	
    	paraMu=0.7
    	paraSigma=15
    	paraC=1
	end-parameters
	

	declarations
		! define the variables for the initialization
		DATA_NUM: integer 	! # of data
	    DATA_NCOL: integer 	! # of columns
	    DATA_N: range 	! range of data
		DATA_COL: range ! range of columns
		_RAWDATA: array(DATA_N,DATA_COL) of real
		
	end-declarations
	
		
	declarations
		
		r_min: real ! minimum value of radius r
		r_max: real ! maximum value of radius r
		r: array(range) of real ! radius array
		
		
		! test variable 
		a: array(range) of real
		
		
		
		! Define the datatype for box
			BOX_CAND=record
				radius: real ! radius
				center: array(DATA_COL) of real ! center point
			end-record
		
			BOX_COL: range
			BOX: array(BOX_COL) of BOX_CAND ! data for boxes
			W: array(DATA_N,BOX_COL) of integer ! incidence matrix
		    
		    BALL: array(1..2) of real ! new ball solution of subproblem ==> BALL(1): pointer of center, BALL(2): radius 
			
		   	! BOX_NUM: integer 
		  	! BOX_N: range ! range of boxes 
	end-declarations
 
    include "common.mos"
    include "subproblem.mos"
    include "result_new_ball.mos"
	
	forward function subproblem(_dualvar: array(range) of real, 
    					_DDATA: array(range,range) of real,
    					_DATA_N: range,_DATA_COL: range,
    					_DATA_NUM: integer,
    					BALL: array(range) of real, ! global variable
    					 _r_max: real, _paraSigma, _paraMu: real) : real
    					 
    					 
    					 
    forward procedure column_gen					 
    					 
    					 
	! Read raw data			
	if(checkfile(D_TYPE)=true) then
			DATA_NUM:=getdatapoints(D_TYPE) 
			DATA_NCOL:=getcolumns(D_TYPE)  
			DATA_N:=1..DATA_NUM    
			DATA_COL:=1..DATA_NCOL    
			getdata(DATA_NCOL,DATA_NUM,D_TYPE,_RAWDATA)
		 
			
	else writeln("There are errors for reading raw data")
		
	end-if	



 	r_min:=calculate_radius_min(_RAWDATA, DATA_NCOL , DATA_NUM)*0.9
 	r_max:=calculate_radius_max(_RAWDATA, DATA_NCOL , DATA_NUM)
	
	
	initialize_box(BOX,_RAWDATA,r_min,DATA_NCOL,DATA_NUM)
	generate_incidence_matrix(W, _RAWDATA, BOX, DATA_NCOL,DATA_N)
    	

   
	


	!-------------------------------------------------------------------------------------
	!  Master Problem
	!-------------------------------------------------------------------------------------
 		declarations
			
			z: array(BOX_COL) of mpvar ! if cube i is selected then zi=1
			sol_z: array(BOX_COL) of real
			xi: array(DATA_N) of mpvar !if xi=0, data k is covered by at least one box. Otherwise, 0
			!kkk: integer ! # of selected boxes
			
			Min_Obj: linctr ! cosntraint for objective function
			Point: array(DATA_N) of linctr ! constraint for point k
			EPS=1e-6
			
			zzz: array(range) of integer ! when output is pointed out
		end-declarations

     
	! set the objective function(density + distance) 
	Min_Obj:= sum(i in BOX_COL)(z(i))+ paraC*sum(k in DATA_N)(xi(k))


		! Set Covering constraint
		forall(k in DATA_N) do
    		Point(k):=sum(i in BOX_COL) (W(k,i)*z(i))  >= 1-xi(k)
		end-do

        ! Constraints for positive data
		forall(i in BOX_COL) z(i) is_binary
	
		setparam('XPRS_VERBOSE', true)
		
		writeln
		writeln

		!Chenck the formulation
		!exportprob("",Obj)

		writeln
		writeln
		
		
		
		
		column_gen
		
		
		
		minimize(Min_Obj)
		
		
		writeln
		
    	writeln
    		writeln("objective value of IP: ", getsol(Min_Obj))
    		writeln
    		
    		_nn:=1
    	 	_kk:=integer(sum(i in BOX_COL) getsol(z(i)))
        	   	 writeln("# of columns: ", getsize(BOX_COL))
        	   	 writeln("# of boxes: ",_kk)
        	   	 writeln
        	   	 writeln("!-------------------------------------------------")
        	   	 
        	 	 
        	forall(i in BOX_COL) do
				if(getsol(z(i))=1) then
			        zzz(_nn):=i
			        _nn+=1	
					writeln("BOX(",i,") is selected")
						
					_kkk:=sum(k in DATA_N) (W(k,i))
					writeln(" # of points in BOX: ",_kkk)
					write("points in the BOX: ")
					
					forall(k in DATA_N) do
						if(W(k,i)=1) then write(k," ")
						end-if
					end-do
					writeln
					writeln("-----------------------------------------------")
					writeln
					
					
				end-if
			end-do	
		


		!write_result_IP( zzz, D_TYPE, paraC)		
	
	
    
	!
    


	!------------------------------------------------------------------------------------------------------------
	
				! Subproblem
		
	!------------------------------------------------------------------------------------------------------------
	    
   function subproblem(_dualvar: array(range) of real, 
    					_DDATA: array(range,range) of real,
    					_DATA_N: range,_DATA_COL: range,
    					_DATA_NUM: integer,
    					BALL: array(range) of real, ! global variable
    					 _r_max: real, _paraSigma, _paraMu: real) : real
		
		declarations
			_sum_dual: real ! objective value of subproblem
			upperbound: real ! upperbound
			lowerbound: real ! lowerbound
			_norm_array: array(range) of real ! temp norm array for each point(center)
			
			_temp: array(range) of real 
			
			
			_temp_ball: array(range, range) of real !i: center pointer, _temp_ball(i,1)=radius, _temp_ball(i,2)=sum of dual variable
			_indcidence_m: array(range) of integer
			
			_N_B: integer ! # of points in a ball
			
		end-declarations
		
		
		
				forall(i in _DATA_N) do
					forall(j in DATA_COL) do 
						_temp(j):=_DDATA(i,j)
					end-do
					 
					max_norm_array(_norm_array,_temp,_DDATA, _DATA_N, DATA_COL)
					
					ii:=1
					while(ii<getsize(_norm_array)+1) do
						lowerbound:=_norm_array(ii)
						get_incidence(_temp, _DDATA, _DATA_N, _indcidence_m, lowerbound)
						_N_B:=sum(k in 1..getsize(_indcidence_m)) _indcidence_m(k)
						
						upperbound:=comparison_min(_r_max,_N_B/_paraSigma)
						
						if(lowerbound<=upperbound) then
							if(centroid_norm(_temp,_RAWDATA,DATA_N,_indcidence_m)<=_paraMu) then
								_temp_ball(i,1):=lowerbound
								_temp_ball(i,2):= sum(k in 1..getsize(_indcidence_m)) _indcidence_m(k)*_dualvar(k)
								ii:=getsize(_norm_array)+1 ! get out of the 'while function'
								
								else ii+=1
									
							end-if
							
							else
								ii+=1
							end-if
							
						
					end-do
				
					if(_temp_ball(i,1)=0) then 
					
						_temp_ball(i,1):=r_min*0.5 ! when the ball should cover the only one point i.e. the center
					end-if
				
				
					
			end-do
		
		
			_sum_dual:=max(i in _DATA_N) _temp_ball(i,2)
		
			forall(i in _DATA_N) do
				if( _temp_ball(i,2)=_sum_dual) then
					BALL(1):=i
					BALL(2):=_temp_ball(i,1)
				end-if
			end-do
		     
		     
		   
			returned:=_sum_dual ! return objective value of subproblem
		

		end-function
	
	!------------------------------------------------------------------------------------------------------------




	!------------------------------------------------------------------------------------------------------------
	
				! CG procedure
		
	!------------------------------------------------------------------------------------------------------------
	 
    procedure column_gen
    	declarations
    		dualpoint: array(DATA_N) of real ! dual variable for each data points
    		new_BALL: array(range) of real ! new_BALL(1): pointer of point & new_BALL(2): radius
    		subobj_1: real ! the value of objective ftn in subproblem -1
    		objval: real ! the value of objective ftn in Master problem
    		
    		bas: basis ! store basis
    		
    		_temp: array(range) of real
    		_w: array(DATA_N) of integer
    		_zzzz: array(range) of integer
    		!_xi: array(range) of integer
    	end-declarations
    	
    	
    	  defcut:=getparam("XPRS_CUTSTRATEGY") ! Save setting of `CUTSTRATEGY' 
 		  setparam("XPRS_CUTSTRATEGY", 0)      ! Disable automatic cuts 
  		  setparam("XPRS_PRESOLVE", 0)         ! Switch presolve off
  		  setparam("zerotol", EPS)             ! Set comparison tolerance of Mosel
 		  nball:=DATA_NUM ! # of columns(ball) ==> initialization of columns
          npass:=1 

 			 while(true) do
   				 minimize(XPRS_LIN,Min_Obj)       ! Solve the LP 

   				 savebasis(bas)                     ! Save the current basis 
    			 objval:= getobjval
    			 
    			 
    			 ! input the solution at each iteration
    			 
    			 forall(i in 1..nball) sol_z(i):=getsol(z(i))   
    			 forall(k in DATA_N) dualpoint(k):=getdual(Point(k))
    			 
    			 subobj_1:=subproblem(dualpoint,_RAWDATA,DATA_N,DATA_COL,DATA_NUM,new_BALL,r_max,paraSigma, paraMu)-1
    			 
    		 
    			 writeln
    			 writeln("-------.", npass,"th iteration-----------")
    			 writeln
    			 
    			 if(subobj_1<=0) then
    			 	writeln("no profitable column found!")
    			 	writeln("objective value: ", objval)
    			 	
    				break
    			 	else
    			 		
    			 		writeln("objective value: ", objval)
    			 		show_new_ball(subobj_1, new_BALL)
    			 		nball+=1
    			 		ii:=integer(new_BALL(1))
    			 		
    			 		! update a new ball information
    			 		BOX(nball).radius:=new_BALL(2)
    			 		forall(j in DATA_COL) do
    			 			BOX(nball).center(j):=_RAWDATA(ii,j)
    			 		end-do
    			 		
    			 		
    			 		create(z(nball))
    			 		z(nball) is_binary
    			 		
    			 		Min_Obj+=z(nball)
    			 		
    			 		
    			 		forall(j in DATA_COL) do 
    			 			_temp(j):=_RAWDATA(ii,j)
    			 		end-do
    			 		
    			 		get_incidence(_temp, _RAWDATA, DATA_N, _w, new_BALL(2))
    			 		
    			 		forall(k in DATA_N) do
    			 			W(k,nball):=_w(k)
    			 		end-do
    			 		forall(k in DATA_N) do
    			 			Point(k)+=_w(k)*z(nball) !add new column
    			 		end-do
    			 		
    			 		loadprob(Min_Obj)
    			 		loadbasis(bas)
    			 	end-if
    			 	
    			 npass+=1		
    	
    		end-do
    	
    	
    	
    		writeln
    		writeln("objective value after column generation: ", objval)
    		writeln
    		
    		nn:=1
    	 	_k:=integer(sum(i in BOX_COL) sol_z(i))
        	   	 writeln("# of columns: ", getsize(BOX_COL))
        	   	 writeln("# of boxes: ",_k)
        	   	 writeln
        	   	 writeln("!-------------------------------------------------")
        	   	 
        	 	 
        	forall(i in BOX_COL) do
				if(sol_z(i)>0) then
			        _zzzz(nn):=i
			        nn+=1	
					writeln("BOX(",i,") is selected")
						
					_kk:=sum(k in DATA_N) (W(k,i))
					writeln(" # of points in BOX: ",_kk)
					write("points in the BOX: ")
					
					forall(k in DATA_N) do
						if(W(k,i)=1) then write(k," ")
						end-if
					end-do
					writeln
					writeln("-----------------------------------------------")
					writeln
					
					
					end-if
			end-do	
        	
        	write_result_LP_relazation( _zzzz, D_TYPE, paraC)	
       	
       	
       	 	
 		
  

  		setparam("XPRS_CUTSTRATEGY", defcut) ! Enable automatic cuts
  		setparam("XPRS_PRESOLVE", 1)         ! Switch presolve on
    	
    	
    	
    	
    	
    end-procedure
    

!------------------------------------------
	!print out new ball
!-----------------------------------------


	
end-model
		
		
		
		
