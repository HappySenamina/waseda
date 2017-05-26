//tansou

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define No_IUnits        4   // input units
#define No_HUnits        3   // hidden units
#define No_OUnits        2   // output units
#define No_Patterns      49   // study patterns
#define S 50

#define Eta              0.4
#define Alpha            0.6
#define ErrorEv          0.0000001 // accepted range of error
#define MaxName          10000

#define fout( x )        (1 / ( 1 + exp( -(x) ) )) 


#define urand()          ( (double)rand() / 10000000000 )

void    forwardPropagation( int p ),backPropagation( int p );
void    testMode(),readFile( FILE  *fp ),initialize();

double out_in[S][No_IUnits],tsignal[S][No_OUnits];
double mi_in[No_IUnits];
double mi_o;
double out_hid[No_HUnits],out_out[No_OUnits];
double witoh[No_HUnits][No_IUnits],dwitoh[No_HUnits][No_IUnits];
double whtoo[No_OUnits][No_HUnits],dwhtoo[No_OUnits][No_HUnits];
double hbias[No_HUnits],dhbias[No_HUnits];
double obias[No_OUnits],dobias[No_OUnits];

int    nlPattern,pTest;

int main(int argc, char* argv[])
{
  int     i=0,j,k;
  double  vError,vErrorM;
  char    sample50[MaxName];
  FILE    *fp;

   fp = fopen("sample50.txt","r");
   if ( fp == NULL )
    { fprintf( stderr,"%s : error!\n",sample50 );
      exit( -1 );
    }
   
  srand( 0 );      
  readFile( fp );  
  initialize();    
  
  vError = ErrorEv + 1.0; 
 
 while(vError > ErrorEv){

      for( j = 0; j < nlPattern; j++ )
        { forwardPropagation( j );
          backPropagation( j );
        }
    
      for( vError = 0.0, j = 0; j < nlPattern; j++ )
        { 
	  
          forwardPropagation( j );

          for( k = 0; k < No_OUnits; k++ )
            vError += pow( tsignal[j][k] - out_out[k], 2.0 );
	    }
      vError *= 0.5;

       if(i%10000==0){
	
	forwardPropagation( nlPattern );
	for( k = 0; k < No_OUnits; k++ )
            vErrorM += pow( tsignal[nlPattern][k] - out_out[k], 2.0 );
	    vErrorM *= 0.5;
	    
	    printf( "%d %.10f %.10f \n",++i,vError,vErrorM );
       }else{
      	i++;
	}
      i++;
       if(i > 10000000){
	 break;
	   }

    
    }

 printf( "\nTest mode ----------------\n" );
 testMode();

  return 0;
}

void forwardPropagation( int p )//forwardPropagation
{
    double sum;

    for(int i = 0; i < No_HUnits; i++ )
      { sum = 0.0;
        for(int j = 0; j < No_IUnits; j++ )
          sum += witoh[i][j]*out_in[p][j];
        out_hid[i] = fout( sum + hbias[i] );
      }
    for(int i = 0; i < No_OUnits; i++ )
      { sum = 0.0;
        for(int j = 0; j < No_HUnits; j++ ){
          sum += whtoo[i][j]*out_hid[j];
        out_out[i] = fout( sum + obias[i] );
	}
      }
}

void backPropagation( int p )//backPropagation
{
    int    i,j;
    double dwih[No_HUnits],dwho[No_OUnits],sum;

    for( i = 0; i < No_OUnits; i++ )
      dwho[i] = ( tsignal[p][i] - out_out[i] )*out_out[i]*( 1.0 - out_out[i] );
    //dwho[i] = ( tsignal[p][i] - out_out[i] )*( 1.0 - out_out[i]*out_out[i] );


    for( i = 0; i < No_HUnits; i++ )
      { for( sum = 0, j = 0; j < No_OUnits; j++ )
          { dwhtoo[j][i] = Eta*dwho[j]*out_hid[i] + Alpha*dwhtoo[j][i];
            whtoo[j][i] += dwhtoo[j][i];
            sum += dwho[j]*whtoo[j][i];
          }
	dwih[i] = out_hid[i]*( 1 - out_hid[i] )*sum;
	//dwih[i] = ( 1 - out_hid[i]*out_hid[i] )*sum;
      }
    for( i = 0; i < No_OUnits; i++ )
      { dobias[i] = Eta*dwho[i] + Alpha*dobias[i];
        obias[i] += dobias[i];
      }


    for( i = 0; i < No_IUnits; i++ )
      { for( j = 0; j < No_HUnits; j++ )
          { dwitoh[j][i] = Eta*dwih[j]*out_in[p][i] + Alpha*dwitoh[j][i];
            witoh[j][i] += dwitoh[j][i];
          }
      }
    for( i = 0; i < No_HUnits; i++ )
      { dhbias[i] = Eta*dwih[i] + Alpha*dhbias[i];
        hbias[i] += dhbias[i];
      }
}

void testMode()
{
  int   i,k;
  char  buf[30];


  for( k = 0; k < nlPattern+1; k++ ){


    rewind( stdin );
    forwardPropagation( k );
            
    /*for( i = 0; i < No_OUnits; i++ ){
      printf( "%d: %f %f \n" ,k+1,out_out[i],tsignal[k][i]);
      }
      printf( "\n" );*/
    printf( "%f %f \n" ,out_out[0],tsignal[k][0]);
  }
printf( "\n" );
  for( k = 0; k < nlPattern+1; k++ ){

    rewind( stdin );
    forwardPropagation( k );
            

    printf( "%f %f \n" ,out_out[1],tsignal[k][1]);
	

  }

}

void readFile( FILE  *fp )
  {
    int   i,j;
    
    fscanf( fp, "%d", &nlPattern );

    for( i = 0; i < nlPattern; i++ )
      { for( j = 0; j < No_IUnits; j++ )
          fscanf( fp, "%lf",&out_in[i][j] );
    for( j = 0; j < No_OUnits; j++ ){
          fscanf( fp, "%lf",&tsignal[i][j] );
    }
      }

    for( i = 0; i < No_IUnits; i++ ){
      fscanf( fp, "%lf",&out_in[nlPattern][i] );
    }for( j = 0; j < No_OUnits; j++ ){
      fscanf( fp, "%lf",&tsignal[nlPattern][j] );
    }

    pTest = nlPattern;
    fclose( fp );
      }
  
void initialize()
  {
    int  i,j;

    for( i = 0; i < No_HUnits; i++ )
      for( j = 0; j < No_IUnits; j++ )
        witoh[i][j] = urand();

    for( i = 0; i < No_OUnits; i++ )
      for( j = 0; j < No_HUnits; j++ )
        whtoo[i][j] = urand();

    for( i = 0; i < No_HUnits; i++ )
      hbias[i] = urand()/10000;
    for( i = 0; i < No_OUnits; i++ )
      obias[i] = urand()/10000;
}
