// Fichero discretelaplacebeltramiC.c
// Recibe una matrices transpuestas
// NOTA: Matlab cuenta desde 1 y C desde 0, por eso a veces se +/- 1

#include "mex.h"
#include <math.h>

#define COLUMNAS_R 1
#define COLUMNAS_P 3
#define EPS 2.220446e-16
#define PI 3.141592653589793


void analizaCeldas(const mxArray *cell_array, int *elemI, mwIndex numCeldas);
void crearSparse(float *lF, double *l, mwIndex *irs, mwIndex *jcs, size_t numnodes);
void discreteLBC(size_t numnodes,int *elemI, float *nodesF, float *eleF, float *kgF, float *avorF, float *lF);
void angleedges(float *u, float *v, float *teta);

//Resta v1-v2 y lo almacena en u
//Resta v3-v2 y lo almacena en v
void restar(float v1[], float v2[],float v3[],float u[],float v[])
{
    int i;
    for (i=0; i<3; i++)
    {
        u[i] = v1[i]-v2[i];
        v[i] = v3[i]-v2[i];
    }
}

//Multiplica v1*v2 y devuelve el resultado
float mult(float v1[], float v2[])
{
    int i;
    float resp = 0;
    for (i=0; i<3; i++)
        resp += v1[i]*v2[i];
    
    return resp;
}

void convert_double2float( double *input_double, float *output_float,size_t size)
{
    int i;
    for (i = 0; i < size; i++)
       output_float[i] = (float) input_double[i];
}

void convert_float2double( float *input_float, double *output_double,size_t size)
{
    int i;
    for (i = 0; i < size; i++)
      output_double[i] = (double) input_float[i];
}

// Esta funcion sera llamada desde Matlab en la forma [Kg,Avor,l]= discretelaplacebeltramiC(ele2node, nodes', elements');
// Todas los inputs matrices se mandan transpuestos
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *nodes, *elements; // inputs; right hand side arg 
	double *kg, *avor, *l; // outputs
    float *nodesF, *eleF, *kgF, *avorF, *lF;
	size_t numnodes, numele, percent_sparse;
    int *elemI;
	mwIndex *irs, *jcs, i;

	// Revisa el numero correcto de argumentos
	if(nrhs != 3) 
		mexErrMsgTxt("Error in number of inputs.");
	else if (nlhs != 3)
		mexErrMsgTxt("Error in number of outputs");

	// Encuentra las dimensiones de los datos
	numnodes = mxGetN(prhs[1]);
	numele = mxGetN(prhs[2]);
	//mexPrintf("numnodes=%d \t numele=%d",numnodes, numele);

	// Crea las matriz para el argumento de retorno
	plhs[0] = mxCreateDoubleMatrix(numnodes, COLUMNAS_R, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(COLUMNAS_R, numnodes, mxREAL);
    percent_sparse = numnodes*6;
    plhs[2] = mxCreateSparse(numnodes,numnodes,percent_sparse,0);

	// Asigna punteros a cada input y output y se hace la conversion de tipos
    elemI = (int*) mxMalloc(percent_sparse*sizeof(int)); //ele2node
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
    
	l = mxGetPr(plhs[2]); 
    irs = mxGetIr(plhs[2]);
    jcs = mxGetJc(plhs[2]);
	lF = (float*) mxMalloc(numnodes*numnodes*sizeof(float));
    
	// Llama la funcion escrita en C
	discreteLBC(numnodes, elemI, nodesF, eleF, kgF, avorF, lF);

    // Hace la conversion de tipos necesaria para los argumentos de retorno
	convert_float2double( kgF, kg, numnodes*COLUMNAS_R);
	convert_float2double( avorF, avor, COLUMNAS_R*numnodes);    
    
	/*for(i=0; i<numnodes;i++)
	{
		mexPrintf("lF[%d] = %f\n", i+1, lF[i] );
	}
*/

    // Crea la matriz sparse l
    crearSparse(lF,l,irs,jcs,numnodes);
    
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
  mwIndex *dim;

  //Cada celda mxArray contiene ( 1 x (5-6) ) mxArray
  for (i=0; i<numCeldas; i++)
  {
      //mexPrintf("\t elementos de la celda: (%d,1) ", i+1);
      cell_element = mxGetCell(cell_array, i);
      if (cell_element != NULL) 
      {
          dim = mxGetDimensions(cell_element);
          col = dim[1];
          //mexPrintf("\t dim(1 X %d) \n", col);
          elem = mxGetPr(cell_element);
          for (j=0; j<col; j++) //dentro de cada celda
          {         
              elemI[i*col+j] =  elem[j]; 
             //mexPrintf("elemF[%d][%d]= %d\n",i+1,j+1, elemI[i*col+j] );
             //mexPrintf("(1,%d) = %f\n", j+1, *elem++ );
          }
      }
      else 
          mexPrintf("celda vacia\n");
      //mexPrintf("-----------------------------------------------\n");
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

//Calcula la curvatura gaussiana KG, Avor y l
void discreteLBC(size_t numnodes,int *elemI, float *nodesF, float *eleF, 
        float *kgF, float *avorF, float *lF)
{
	int i,h,b,j,k,dimension,temp2;
	float bf1, bf2, bf3, alpha, beta, thetat, temp;
	float vi[3], vj[3], vk[3], u[3], v[3];
    float *dist2;
    
    dist2 = (float*) mxMalloc(numnodes*sizeof(float));    
	for(i=0; i<numnodes; i++)
	{
		thetat = 0;	
        //memset(dist2,0,numnodes*sizeof(float)); //sin esto es MUCHO + lento
        temp2 = numnodes*i;
		if(i>11)
            dimension=6;
        else
            dimension=5;
		for(h=0;h<dimension;h++) 
		{
			//b is elements adyacent to i node
			b = elemI[i*dimension+h] - 1;  
            b=b*3;
			bf1 = eleF[b+0];
			bf2 = eleF[b+1];
			bf3 = eleF[b+2];
			//mexPrintf(" b=%d ______ bf1=%f, bf2=%f, bf3=%f",b+1,bf1,bf2,bf3);

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
            b=i*3;
			vi[0] = nodesF[b+0];
			vi[1] = nodesF[b+1];
			vi[2] = nodesF[b+2];
            b=j*3;
			vj[0] = nodesF[b+0];
			vj[1] = nodesF[b+1];
			vj[2] = nodesF[b+2];
            b=k*3;
			vk[0] = nodesF[b+0];
			vk[1] = nodesF[b+1];
		 	vk[2] = nodesF[b+2];

			//mexPrintf("\n vi[0]=%f,vi[1]=%f, vi[2]=%f",vi[0],vi[1],vi[2]);
			//mexPrintf("\n vj[0]=%f,vj[1]=%f, vj[2]=%f",vj[0],vj[1],vj[2]);
			//mexPrintf("\n vk[0]=%f,vk[1]=%f, vk[2]=%f\n",vk[0],vk[1],vk[2]);
					
			//angles
            restar(vi,vk,vj,u,v);			        
			angleedges(u, v, &beta); // (vi-vk,vj-vk)
            restar(vi,vj,vk,u,v); 
			angleedges(u, v, &alpha); // (vi-vj,vk-vj)
            
            temp=0;
            for (b=0; b<3; b++)
                temp+= pow(u[b],2);
            dist2[j] = temp/8;
                    
            restar(vk,vi,vj,u,v);
            angleedges(u, v, &temp); // (vk-vi,vj-vi)
			thetat = thetat + temp;
			//mexPrintf("\nbeta=%f,alpha=%f,thetat=%f\n",beta,alpha,thetat);
			
            //add weight
            temp = 1/tan(beta);
            lF[temp2+j]+= temp;
            temp = 1/tan(alpha);
            lF[temp2+k]+= temp;
            //mexPrintf("------------------[%d][%d]---------------------\n",i+1, h+1);
		}
        temp=0;
        for(b=0; b<numnodes; b++)
        {
            if( dist2[b]!= 0)
            {
                temp += lF[temp2+b]*dist2[b];
                dist2[b] = 0;
            }
        }        
    	avorF[i] = temp;
        
        //Gaussian Curvature
		kgF[i] = (2*PI - thetat)/temp;
	}

	mxFree(dist2);
}

//Encuentra el angulo teta
void angleedges(float *u, float *v, float *teta)
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
	*teta = (float) acos(mult(u,v)/(du*dv));
}
