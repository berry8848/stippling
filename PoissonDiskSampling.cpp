#include "PoissonDiskSampling.h"



PoissonDiskSampling::PoissonDiskSampling(
    int         width,
    int         height,
    double      pitch
)
{
    width_  = width;
    height_ = height;
    pitch_  = pitch;

    // �m�[�h�z��̗p��
    nodes.resize( width_ * height_ );

    // �����m�[�h�̐ݒ�
    for ( int j=0; j<height_; j++ ) {
        for( int i=0; i<width_; i++ ) {        
            nodes[i+j*width] = Node( i, j, 1.0 );
        }
    }

    // ���������֐��̏�����
    srand( (unsigned int)( time(0) ) );

    f = Converter();
}



PoissonDiskSampling::~PoissonDiskSampling(){

}



void PoissonDiskSampling::sample(
    std::vector<Eigen::Vector3d>&   points
)
{
    std::vector<Node>::iterator n_it;
    std::vector<Node>::iterator n_itN;
    std::vector<int> neighbors; 
    std::deque<Eigen::Vector3d> active_list;

    // �_��z��̏�����
    points.resize( 0 );
 
    // �����_�̐���
    double xx = (double)width_  - 2.0;
    double yy = (double)height_ - 2.0;
    double zz = 0.0;
    Eigen::Vector3d init_point = Eigen::Vector3d( xx, yy, zz );

    // ���_���X�g�ɒǉ�
    active_list.push_back( init_point );

    // �m��_�Ƃ��ĕێ�
    insertPoint( nodes, points, init_point );

    // ���_�����O�ɂȂ�܂ŏ������J��Ԃ�
    while( active_list.size() != 0 ) {

        // ���_
        Eigen::Vector3d center = active_list[0];
        bool is_valid = false;

        // ���_��������m�[�h�̎擾
        int pid_x = (int)( center[0] / pitch_ );
        int pid_y = (int)( center[1] / pitch_ );
        int n_id  = pid_x + pid_y * width_;
        n_it = nodes.begin() + n_id;

        // �m�[�h���ێ�����W������_�ԋ���������
        double rr = f( n_it->val );        

        // rr �` rr * 2.0 �͈̔͂ɐV�����_��ǉ�
        int add_npnt = 10;

        for ( int ipnt=0; ipnt<add_npnt; ipnt++ ) {

            // �͈͓��̒ǉ��_�̈ʒu�������_���Ɍ���
            double  lower_rr = rr;
            double  upper_rr = rr * 2.0;
            Eigen::Vector3d add_pnt = genRandomPoint( center, lower_rr, upper_rr );

            // �ǉ��_��������m�[�h�̎擾
            int pid_xN = (int)( add_pnt[0] / pitch_ );
            int pid_yN = (int)( add_pnt[1] / pitch_ );
            int n_idN  = pid_xN + pid_yN * width_;
            n_itN =nodes.begin()+ n_idN;

            // �m�[�h���ێ�����W������_�ԋ���������
            double rrN = f( n_itN->val );

            // �ǉ��_�̓_�ԋ������ɑ��̓_�����݂��Ȃ����m�F
            bool is_exist = existNeighbors( add_pnt, rrN, points );
            if ( is_exist == false ){

                // ���_���X�g�ɒǉ�
                active_list.push_back( add_pnt );

                // �m��_�Ƃ��ĕێ�
                insertPoint( nodes, points, add_pnt );

                is_valid = true;
            }
        }       

        // �ǉ��_���Ȃ��ꍇ�͒��ړ_�����_����O��
        if ( is_valid == false ) {      
            active_list.pop_front();
        }
    }
}



void PoissonDiskSampling::insertPoint(
    std::vector<Node>&              nodes,
    std::vector<Eigen::Vector3d>&   points,
    const Eigen::Vector3d&          point
)
{
    // �m��_��������m�[�h�C���f�b�N�X�̎Z�o
    int pid_x = (int)( point[0] / pitch_ );
    int pid_y = (int)( point[1] / pitch_ );
    int n_id  = pid_x + pid_y * width_;

    // �m��_���������m�[�h�Ƀ��x��
    nodes.at( n_id ).visited = true;
    
    // �m��_�̕ێ�
    points.push_back( point );

    // �m��_�̔z����C���f�b�N�X���m�[�h�ɕێ�
    int index = (int)( points.size() - 1 );
    nodes.at(n_id).p_id.push_back( index );
}



void PoissonDiskSampling::setDensityFunc(
    const std::vector<double>&  dfunc
)
{
    // �C���^�[�i���`�F�b�N
    assert( (int)dfunc.size() == width_*height_ );

    // �e�m�[�h�ɌW�����Z�b�g
    for ( int ii=0; ii<nodes.size(); ii++ ) {
        nodes[ii].val = dfunc[ii];
    }
}



void PoissonDiskSampling::setConverter(boost::function<double(double)> _f){ 
    f = _f;
}



Eigen::Vector3d PoissonDiskSampling::genRandomPoint(
    const Eigen::Vector3d&  point,
    double                  low_rr,
    double                  upp_rr
)
{
    Eigen::Vector3d     p;
    
    double low_rr2 = low_rr * low_rr;
    double upp_rr2 = upp_rr * upp_rr; 

    while ( 1 ) {
        // point �����_�Ƃ������W�l�Ń����_���ɓ_�𔭐�
        double      fittol;
        double      dx,dy,xcoef,ycoef,xrand,yrand,max_rand;
        
        max_rand = (double)RAND_MAX;
        fittol   = max_rand * 0.001;

        xrand = (double)rand();  // 0 - RAND_MAX
        xcoef = ( xrand * 2.0 - max_rand + fittol ) / ( max_rand + 2.0 * fittol );
        dx    = xcoef * upp_rr;

        yrand = (double)rand();  // 0 - RAND_MAX
        ycoef = ( yrand * 2.0 - max_rand + fittol ) / ( max_rand + 2.0 * fittol );
        dy    = ycoef * upp_rr;

        //double dx = ((double)rand()+1.0)/((double)RAND_MAX+2.0)*2*upp_rr - upp_rr;
        //double dy = ((double)rand()+1.0)/((double)RAND_MAX+2.0)*2*upp_rr - upp_rr;
        //double dz = ((double)rand()+1.0)/((double)RAND_MAX+2.0)*2*upp_rr - upp_rr;  

        // �~�͈̔͊O�̓X�L�b�v���čăg���C
        double cur_rr2 = dx * dx + dy * dy;
        if ( cur_rr2 < low_rr2 || upp_rr2 < cur_rr2 ) {
            continue;
        }

        double x = point[0] + dx;
        double y = point[1] + dy;

        // �L���͈͊O�̓X�L�b�v���čăg���C
        if ( x < 0  ||  width_  < x )  continue;
        if ( y < 0  ||  height_ < y )  continue;

        // �ʒu�m��
        p = Eigen::Vector3d( x, y, 0.0 );
        break;
    }

    return p;
}



bool PoissonDiskSampling::existNeighbors(
    Eigen::Vector3d&    center,
    double              rr,
    std::vector<Eigen::Vector3d>&   points
)
{
    bool exist = false;

    double dx  = center[0];
    double dy  = center[1];
    double rr2 = rr * rr;

    // ���ړ_���甼�a�͈͓��̃m�[�h�̒T��
    std::vector<int> near_nodes;

    int min_x = (int)( ( dx - rr ) / pitch_ );
    int max_x = (int)( ( dx + rr ) / pitch_ );
    int min_y = (int)( ( dy - rr ) / pitch_ );
    int max_y = (int)( ( dy + rr ) / pitch_ );

    for ( int j=min_y; j<=max_y; j++ ) {
        if ( j <  0      ) continue;
        if ( j >= height_ ) continue;

        for ( int i=min_x; i<=max_x; i++ ) {
            if ( i <  0      ) continue;
            if ( i >= width_ ) continue;

            int n_id  = i + j * width_;
            near_nodes.push_back( n_id );
        }
    }

    // �d���m�[�h�̍폜
    int         num_node,num_remove,ii,jj;
    int         ii_id,jj_id;

    num_node   = (int)( near_nodes.size() );
    num_remove = 0;

    for ( ii=0; ii<num_node-1; ii++ ) {
        ii_id = near_nodes[ii];
        if ( ii_id == -1 )  continue;

        for ( jj=ii+1; jj<num_node; jj++ ) {
            jj_id = near_nodes[jj];
            if ( jj_id == -1 )  continue;

            if ( ii_id == jj_id ) {
                num_remove ++;
                near_nodes[jj] = -1;
            }
        }
    }

    // �m�[�h�h�c�̃\�[�g
    for ( ii=0; ii<num_node-1; ii++ ) {
        for ( jj=ii+1; jj<num_node; jj++ ) {
            ii_id = near_nodes[ii];
            jj_id = near_nodes[jj];

            if ( jj_id > ii_id ) {
                 near_nodes[ii] = jj_id;
                 near_nodes[jj] = ii_id;
            }
        }
    }

    num_node = num_node - num_remove;

    // �m�[�h���̓_�Q�Ŕ��a�����߂Â��_���Ȃ����m�F
    bool        visited;
    int         num_pnts,ipnt,nodidx,pntidx;
    double      dist2;
    Eigen::Vector3d  point;

    for ( ii=0; ii<num_node; ii++ ) {
        nodidx = near_nodes[ii];

        visited = nodes[nodidx].visited;
        if ( visited == false )  continue;

        num_pnts = (int)( nodes[nodidx].p_id.size() );

        for ( ipnt=0; ipnt<num_pnts; ipnt++ ) {
            pntidx = nodes[nodidx].p_id[ipnt];
            point  = points[pntidx];

            dx = point[0] - center[0];
            dy = point[1] - center[1];

            dist2 = dx * dx + dy * dy;
            if ( dist2 < rr2 ) {
                exist = true;
                break;              
            }
        }

        if ( exist == true )  break;
    }

    return exist;
}
