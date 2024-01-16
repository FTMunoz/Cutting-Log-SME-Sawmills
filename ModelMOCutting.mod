/*********************************************
Model: IBM ILOG CPLEX Optimization Studio
 *********************************************/

 // Lectura  de Datos
int diametros = ...;   
int patrones = ...;     
int productos = ...; 

range r_diametros = 1..diametros;
range r_patrones = 1..patrones;
range r_productos= 1..productos;


// Lectura Hoja Matriz
string mat_diametro[1..patrones] = ...;
string mat_patron[1..patrones] = ...;
int mat_cantprod[1..patrones][1..productos] =...;
float    mat_perd[1..patrones] = ...;

 // Lectura Hoja Dda
string dem_prod[1..productos] = ...;
int    dem_cant[1..productos] = ...;
float    dem_w[1..productos] = ...;

 // Lectura Hoja Disp
string disp_diametro[1..diametros] = ...;
int    disp_cantidad[1..diametros] = ...;

// Variables de Decisión
dvar int+ x[1..patrones];
dvar float+ n[1..productos];
dvar float+ FO1;
dvar float+ FO2;

minimize  FO1+FO2; 
												
subject to {	

	FO1 == sum(i in 1..patrones) mat_perd[i]*x[i];
	
	FO2 == sum(p in 1..productos) (n[p]*dem_w[p]) ;
	
	forall(d in 1..diametros) 
		DISPONIBILIDAD:
		sum(i in 1..patrones: mat_diametro[i]==disp_diametro[d]) x[i] <= disp_cantidad[d] ;

	forall(p in 1..productos) 
		DEMANDA:
		sum(i in 1..patrones) mat_cantprod[i][p]*x[i] >= dem_cant[p] ;	

	forall(p in 1..productos) 
		DESV:
		n[p] == sum(i in 1..patrones) mat_cantprod[i][p]*x[i] - dem_cant[p] ;	
		
	
}


 
 
main									
{	cplex.epgap = 0.01;
	var ptos=25;	
	var epsilon=0;
	
	var R = new IloOplOutputFile("Resumen_Pareto.txt");	
	R.writeln("Objetivo1\tObjetivo2");	
	var S = new IloOplOutputFile("Salida.txt");	
	S.writeln("Objetivo1\tObjetivo2\t(patron;diametro;cantidad)....");
	
   thisOplModel.generate();
   //thisOplModel.convertAllIntVars();
   
   thisOplModel.FO1.LB=0;
   thisOplModel.FO1.UB=Infinity;
   thisOplModel.FO2.LB=0;
   thisOplModel.FO2.UB=Infinity;
   cplex.setObjCoef(thisOplModel.FO1,1);
   cplex.setObjCoef(thisOplModel.FO2,1e-12);
   cplex.solve();
   //R.writeln(thisOplModel.FO1,"\t",thisOplModel.FO2);
   writeln(thisOplModel.FO1,"\t",thisOplModel.FO2);
   S.writeln("");
   S.write(thisOplModel.FO1,"\t",thisOplModel.FO2,"\t");
   for (var i in thisOplModel.r_patrones)							
	{	if(thisOplModel.x[i]!=0) S.write("(", thisOplModel.mat_patron[i], ";", thisOplModel.mat_diametro[i], ";",thisOplModel.x[i],") ");
   	} 
   var f1min = thisOplModel.FO1.solutionValue; 
   var f2max = thisOplModel.FO2.solutionValue; 
   writeln("Modelo Obj1 OK");	
   
   thisOplModel.FO1.LB=0;
   thisOplModel.FO1.UB=Infinity;
   thisOplModel.FO2.LB=0;
   thisOplModel.FO2.UB=Infinity;
   cplex.setObjCoef(thisOplModel.FO1,1e-12); // muy chico
   cplex.setObjCoef(thisOplModel.FO2,1);
   cplex.solve();
   R.writeln(thisOplModel.FO1,"\t",thisOplModel.FO2);
   writeln(thisOplModel.FO1,"\t",thisOplModel.FO2);
   S.writeln("");
   S.write(thisOplModel.FO1,"\t",thisOplModel.FO2,"\t");
   for (var i in thisOplModel.r_patrones)							
	{	if(thisOplModel.x[i]!=0) S.write("(", thisOplModel.mat_patron[i], ";", thisOplModel.mat_diametro[i], ";",thisOplModel.x[i],") ");
    } 
   var f2min = thisOplModel.FO2.solutionValue;
   var f1max = thisOplModel.FO1.solutionValue;
   
   epsilon = (f2max - f2min)/ptos;
   
   writeln("Modelo Obj2 OK");
   writeln("Frontera de pareto (Obj1, Obj2)");
   
   var REV = new IloOplOutputFile("Rev.txt");	
   

   var k=1;
   cplex.setObjCoef(thisOplModel.FO1,1);
   cplex.setObjCoef(thisOplModel.FO2,0);
   while(k<ptos)
   {	
   		//thisOplModel.FO2.LB=0;
   		thisOplModel.FO2.UB = f2min + epsilon*k;
   		
   		cplex.solve();
   		//obj1 = thisOplModel.FO1;
   		//beta = beta + epsilon;
   		//antes de imprimir preguntar si es el mismo punto
   			R.writeln(thisOplModel.FO1,"\t",thisOplModel.FO2);
	   		writeln(thisOplModel.FO1,"\t",thisOplModel.FO2);
	   		S.writeln("");
	   		S.write(thisOplModel.FO1,"\t",thisOplModel.FO2,"\t");
	   		for (var i in thisOplModel.r_patrones)							
			{	if(thisOplModel.x[i]!=0) S.write("(", thisOplModel.mat_patron[i], ";", thisOplModel.mat_diametro[i], ";",thisOplModel.x[i],") ");
	       	} 
	       	
	       	REV.writeln(thisOplModel.FO1,"\t",thisOplModel.FO2);
	       	for (var i in thisOplModel.r_patrones)							
			{	REV.write(thisOplModel.x[i]," ");
	       	} 
	       	REV.writeln("");
	       	for (var i in thisOplModel.r_productos)							
			{	REV.write(thisOplModel.n[i]," ");
	       	} 
	       	REV.writeln("");
       
       	k=k+1;
   }      
    
//    	thisOplModel.FO1.LB=0;
//   		thisOplModel.FO1.UB=alfa;
//   		thisOplModel.FO2.LB=0;
//   		thisOplModel.FO2.UB=Infinity;
//   		cplex.setObjCoef(thisOplModel.FO1,0);
//  		cplex.setObjCoef(thisOplModel.FO2,1);
//   		cplex.solve();
//   		R.writeln(thisOplModel.FO1,"\t",thisOplModel.FO2);
//   		writeln(thisOplModel.FO1,"\t",thisOplModel.FO2);
//   		S.writeln("");
//   		S.write(thisOplModel.FO1,"\t",thisOplModel.FO2,"\t");
//   		for (var i in thisOplModel.r_patrones)							
//		{	if(thisOplModel.x[i]!=0) S.write("(", thisOplModel.mat_patron[i], ";", thisOplModel.mat_diametro[i], ";",thisOplModel.x[i],") ");
//       	} 
    R.writeln(f1min,"\t",f2max);
	
   	
	R.close();
	S.close();
	REV.close();
						
}																						