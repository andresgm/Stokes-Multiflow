// Fichero discreteLBC.cpp
// Recibe una matrices transpuestas
// NOTA: Matlab cuenta desde 1 y C desde 0, por eso a veces se +/- 1

#include "mex.h"
#include <math.h>
#include <windows.h>

#define COLUMNAS_R 1
#define COLUMNAS_P 3
#define EPS 2.220446e-16
#define PI 3.141592653589793


void analizaCeldas(const mxArray *cell_array, int *elemI, mwIndex numCeldas);
void crearSparse(float *lF, double *l, mwIndex *irs, mwIndex *jcs, size_t numnodes);
void discreteLBC(size_t numnodes,int *elemI, float *nodesF, float *eleF, float *kgF, float *avorF, float *lF);
void angleedgesCPP(float *u, float *v, float *teta);

//Resta v1-v2 y lo almacena en u
//Resta v3-v2 y lo almacena en v
void restarCPP(float v1[], float v2[],float v3[],float u[],float v[])
{
    int i;
    for (i=0; i<3; i++)
    {
        u[i] = v1[i]-v2[i];
        v[i] = v3[i]-v2[i];
    }
}

//Multiplica v1*v2 y devuelve el resultado
float multCPP(float v1[], float v2[])
{
    int i;
    float resp = 0;
    for (i=0; i<3; i++)
        resp += v1[i]*v2[i];
    
    return resp;
}

//Transforma los valores de doble precision a precision simple
void convert_double2float( double *input_double, float *output_float,size_t size)
{
    int i;
    for (i = 0; i < size; i++)
       output_float[i] = (float) input_double[i];
}

//Transforma los valores de simple precision a precision doble
void convert_float2double( float *input_float, double *output_double,size_t size)
{
    int i;
    for (i = 0; i < size; i++)
      output_double[i] = (double) input_float[i];
}

// Esta funcion sera llamada desde Matlab en la forma [Kg,Avor,L]= discreteLBC(ele2node, nodes', elements')
// Todas los inputs matrices se mandan transpuestos
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *nodes, *elements; // inputs; right hand side arg 
	double *kg, *avor, *L; // outputs
    float *nodesF, *eleF, *kgF, *avorF, *lF;
	size_t numnodes, numele;
    int *elemI;
	mwIndex *irs, *jcs, i; //valores para la SparceMatrix
	
	// Revisa el numero correcto de argumentos
	if(nrhs != 3) 
		mexErrMsgTxt("Error in number of inputs.");
	else if (nlhs != 3)
		mexErrMsgTxt("Error in number of outputs");

	// Encuentra las dimensiones de los datos
	numnodes = mxGetN(prhs[1]);
	numele = mxGetN(prhs[2]);

	// Crea las matriz para el argumento de retorno
	plhs[0] = mxCreateDoubleMatrix(numnodes, COLUMNAS_R, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(COLUMNAS_R, numnodes, mxREAL);
    plhs[2] = mxCreateSparse(numnodes,numnodes,numnodes*7,mxREAL);

	// Asigna punteros a cada input y output y se hace la conversion de tipos
    elemI = (int*) mxMalloc(numnodes*6*sizeof(int)); //ele2node
    analizaCeldas(prhs[0], elemI, numnodes);     
	nodes = mxGetPr(prhs[1]);
	nodesF = (float*) mxMalloc(numnodes*COLUMNAS_P*sizeof(float));
	convert_double2float( nodes, nodesF, numnodes*COLUMNAS_P);
	elements = mxGetPr(prhs[2]);
	eleF = (float*) mxMalloc(numele*COLUMNAS_P*sizeof(float));	
	convert_double2float( elements, eleF, numele*COLUMNAS_P);

	kg = mxGetPr(plhs[0]);
	kgF = (float*) mxMalloc(numnodes*COLUMNAS_R*sizeof(float));
	avor = mxGetPr(plhs[1]);
	avorF = (float*) mxMalloc(COLUMNAS_R*numnodes*sizeof(float));    
	L = mxGetPr(plhs[2]); 
    irs = mxGetIr(plhs[2]);
    jcs = mxGetJc(plhs[2]);
	lF = (float*) mxMalloc(numnodes*numnodes*sizeof(float));

	// Llama la funcion escrita en C
	discreteLBC(numnodes, elemI, nodesF, eleF, kgF, avorF, lF);
 
    // Hace la conversion de tipos necesaria para los argumentos de retorno
	convert_float2double( kgF, kg, numnodes*COLUMNAS_R);
	convert_float2double( avorF, avor, COLUMNAS_R*numnodes);   

    // Crea la matriz sparse L
    crearSparse(lF,L,irs,jcs,numnodes);
    
    // Libera el espacio reservado
	mxFree(nodesF);
	mxFree(eleF);
	mxFree(kgF);
	mxFree(avorF); 
	mxFree(lF);   
	mxFree(elemI);
}

//Extrae los datos del arreglo de celdas (numceldas x 1)
void analizaCeldas(const mxArray *cell_array, int *elemI, mwIndex numCeldas)
{
  const mxArray *cell_element;
  double *elem;
  mwIndex i, j, col;
  //Cada celda mxArray contiene ( 1 x (5-6) ) mxArray
  for (i=0; i<numCeldas; i++)
  {
      cell_element = mxGetCell(cell_array, i);
      if (cell_element != NULL) 
      {
          col = mxGetNumberOfElements(cell_element);
          elem = mxGetPr(cell_element);
          for (j=0; j<col; j++) //dentro de cada celda
              elemI[i*col+j] =  elem[j];    
      }
      else 
          mexPrintf("celda vacia\n");
  }
}

//Crea una matriz sparse desde una normal
//irs apunta a un arreglo de indices de la fila
//jcs apunta a un arreglo de indices de la columna
void crearSparse(float *lF, double *l, mwIndex *irs, mwIndex *jcs, size_t numnodes)
{
    mwSize i;
	mwIndex	j,k;
    k = 0; 
    
    //Copy nonzeros
    for (j=0; j<numnodes; j++) 
	{
      jcs[j] = k;
      for (i=0; i<numnodes ; i++)
	  {
		if (lF[i] != 0.0) 
		{
			l[k] = (double) lF[i];
			irs[k] = i;
			k++;
        }    
      }
      lF += numnodes;
    }
    jcs[numnodes] = k;
}

//Calcula la curvatura gaussiana KG, Avor y L
void discreteLBC(size_t numnodes,int *elemI, float *nodesF, float *eleF, 
        float *kgF, float *avorF, float *lF)
{
	int i,h,b,j,k,dimension,temp2;
	float bf1, bf2, bf3, alpha, beta, thetat, temp, temp3;
	float vi[3], vj[3], vk[3], u[3], v[3];
    float *dist2;
    
    dist2 = (float*) mxMalloc(numnodes*sizeof(float));    
	for(i=0; i<numnodes; i++)
	{
		thetat = 0;	
        memset(dist2,0,numnodes*sizeof(float)); //necesario para accelerar el calculo
        temp2 = numnodes*i;
		if(i>11)
            dimension=6;
        else
            dimension=5;
		for(h=0;h<dimension;h++) 
		{
			//b is elements adyacent to i node
			b = elemI[i*dimension+h] - 1;  
			bf1 = eleF[b*3+0];
			bf2 = eleF[b*3+1];
			bf3 = eleF[b*3+2];

			//define edges
			if(bf1 == i+1)
			{
				j = bf2-1;
				k = bf3-1;
			}
			else if(bf2 == i+1)
			{
				j = bf3-1;
				k = bf1-1;
			}
			else //(bf3 == i+1)
			{
				j = bf1-1;
				k = bf2-1;
			}
            
			vi[0] = nodesF[i*3+0];
			vi[1] = nodesF[i*3+1];
			vi[2] = nodesF[i*3+2];

			vj[0] = nodesF[j*3+0];
			vj[1] = nodesF[j*3+1];
			vj[2] = nodesF[j*3+2];

			vk[0] = nodesF[k*3+0];
			vk[1] = nodesF[k*3+1];
		 	vk[2] = nodesF[k*3+2];
					
			//angles
            restarCPP(vi,vk,vj,u,v);			        
			angleedgesCPP(u, v, &beta); // (vi-vk,vj-vk)
            restarCPP(vi,vj,vk,u,v); 
			angleedgesCPP(u, v, &alpha); // (vi-vj,vk-vj)
            
            temp=0;
            for (b=0; b<3; b++)
                temp+= pow(u[b],2);
            dist2[j] = temp/8;
                    
            restarCPP(vk,vi,vj,u,v);
            angleedgesCPP(u, v, &temp); // (vk-vi,vj-vi)
			thetat = thetat + temp;
			
            //add weight
            temp = 1/tan(beta);
            lF[temp2+j]+= temp;
            temp = 1/tan(alpha);
            lF[temp2+k]+= temp;
		}
        temp=0;
        temp3=0;
        for(b=0; b<numnodes; b++)
        {
            temp3 += lF[temp2+b];
            temp += lF[temp2+b]*dist2[b];
        }        
    	avorF[i] = temp;
        lF[temp2+i] = -temp3;
        //Gaussian Curvature
		kgF[i] = (2*PI - thetat)/temp;
	}

	mxFree(dist2);
}

//Encuentra el angulo teta
void angleedgesCPP(float *u, float *v, float *teta)
{
    float du = 0;
    float dv = 0;
    int i;
    
    for (i=0; i<3; i++)
    {
        du += pow(u[i],2);
        dv += pow(v[i],2);
    }
    du = sqrt(du);
    dv = sqrt(dv);
    if (du < EPS)
        du = (float) EPS;
    if (dv < EPS)
        dv = (float) EPS;
    
    //cos(teta) = <u,v>/(|u|*|b|)
	*teta = (float) acos(multCPP(u,v)/(du*dv));
 }