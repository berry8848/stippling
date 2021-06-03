#include    <iostream>
#include    <string>
#include    <fstream>
#include    <opencv2/opencv.hpp>
#include    <Eigen/Core>
#include    "PoissonDiskSampling.h"
//#include "GeomRenderer.h"
#include    "IO.h"



// l = √3d((1+√(1-σ))/σ)
class F : public PoissonDiskSampling::Converter{
    public: 
        double operator()( double val ) {
            if ( fabs(val) < thresh() ) {
                return max();
            }
            return  sqrt( 3.0 ) * 0.1 * ( 1.0 + sqrt( 1.0 - val ) ) / val;
        };
        static double thresh(){
            return 0.0001;
        }
        static double max(){
            return 1000;
        }
    };


// 説  明    ：入力した画像を点群に換算
// 入力引数１：画像ファイル名
// 入力引数２：点群ファイル名  

int main( int argc, char* argv[] )
{
    /*----- Declare and Initialize -----------------------------------*/
    long            rc  = 0;            /* return code                */
    long            frc = 0;            /* function return code       */

    char            errmsg[128];
    char            img_fname[512],pnt_fname[512];
    long            num_pnts,ipnt,img_ww,img_hh,ii,jj,size,index;
    double          val,min_val,max_val,ww_max,ww_min,coef;

    Eigen::Vector3d pnt;
    std::vector<Eigen::Vector3d> points;


    FILE*           fp_err = NULL;
    FILE*           fp_ply = NULL;


    /*----- Check Input Argument -----*/
    if ( argc != 3 ) {
        rc = 1;
        goto ERR;
    }


    /*----- Get File Name --------------------------------------------*/
    strcpy( img_fname, argv[1] );
    strcpy( pnt_fname, argv[2] );


    {
        /*----- Read Image File --------------------------------------*/
        cv::Mat im = cv::imread( img_fname, 0 );

        img_ww = im.cols;
        img_hh = im.rows;

        std::vector<double>  ww_vec( img_ww * img_hh );

        for ( ii=0; ii<img_hh; ii++ ) {
            for ( jj=0; jj<img_ww; jj++ ) {     
                 val = (double)im.at<unsigned char>( ii, jj );
                 ww_vec[ ii * img_ww + jj ] = val;
                 //w_vec[j*width+i] = (val)/(double)255*(w_max-w_min)+w_min;
            }
        }

        size = img_ww * img_hh;

        for ( ii=0; ii<size; ii++ ) {
            val = ww_vec[ii];

            if ( ii == 0 ) {
                min_val = val;
                max_val = val;
            }
            else {
                if ( val < min_val )  min_val = val;
                if ( val > max_val )  max_val = val;
            }
        }

        ww_min = 0.5;
        ww_max = 5.0;

        for ( ii=0; ii<img_hh; ii++ ) {
            for ( jj=0; jj<img_ww; jj++ ) {     
                 index = ii * img_ww + jj;
                 val   = ww_vec[index];
                 coef  = ( val - min_val ) / ( max_val - min_val );
                 ww_vec[index] = ( ww_max - ww_min ) * coef + ww_min;
            }
        }


        /*----- Poisson Disk Sampling --------------------------------*/
        PoissonDiskSampling pds( img_ww, img_hh, 1.0 );
        pds.setDensityFunc( ww_vec );
        pds.sample( points );
    }

    num_pnts = (long)( points.size() );
    if ( num_pnts < 1 ) {
        rc = 2;
        goto ERR;
    }


    /*----- Write Points to File ( PLY Format ) ----------------------*/
    fp_ply = fopen( "POINTS.ply", "w" );
    if ( fp_ply == NULL ) {
        rc = 3;
        goto ERR;
    }

    fprintf( fp_ply, "ply \n" );
    fprintf( fp_ply, "format ascii 1.0 \n" );
    fprintf( fp_ply, "element vertex %d \n", num_pnts );
    fprintf( fp_ply, "property double x \n" );
    fprintf( fp_ply, "property double y \n" );
    fprintf( fp_ply, "property double z \n" );
    fprintf( fp_ply, "end_header \n" );
    
    for ( ipnt=0; ipnt<num_pnts; ipnt++ ) {
        pnt = points[ipnt];
        fprintf( fp_ply, "%lf %lf %lf \n", pnt[0], pnt[1], pnt[2] );
    }

    fclose( fp_ply );
    fp_ply = NULL;


    goto OWARI;
    /*----- Error Message --------------------------------------------*/
ERR:
    switch( rc ){
        case 1 :
            sprintf( errmsg, "no. of arguments not 3. " );
            break;
        case 2 :
            sprintf( errmsg, "no. of point < 1. " );
            break;
        case 3 :
            sprintf( errmsg, "err in fopen(). " );
            break;
        default:
            sprintf( errmsg, "unknown error" );
            break;
    }

    fp_err = fopen( "ERRLOG.txt", "w" );
    if ( fp_err != NULL ) {
        fprintf( fp_err, "%s", errmsg );
        fclose( fp_err );
        fp_err = NULL;
    }


    /*----- Closing and Return ---------------------------------------*/
OWARI:
    /*----- Free -----*/
    if ( fp_err != NULL )  fclose( fp_err );
    if ( fp_ply != NULL )  fclose( fp_ply );

    //printf("\n");
    //printf( ">>>>> Hit Any Key + RET to Continue ! " );
    //scanf( "%c", &cdummy );

    /*----- Terminate -----*/
    return( rc );
}