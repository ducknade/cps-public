// vim: set ts=2 sw=2 expandtab:
#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  LinkBuffer class methods and some Lattice class methods involving
  the link buffer.

  $Id $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008/09/18 15:23:17 $
//  $Header: /space/cvs/cps/cps++/src/util/lattice/lattice_base/link_buffer.C,v 1.8 2008/09/18 15:23:17 chulwoo Exp $
//  $Id: link_buffer.C,v 1.8 2008/09/18 15:23:17 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $Revision: 1.8 $
//  $Source: /space/cvs/cps/cps++/src/util/lattice/lattice_base/link_buffer.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include <util/vector.h>
#include <util/verbose.h>
#include <util/time_cps.h>
#include <comms/nga_reg.h>
#include <comms/scu.h>
#include <comms/cbuf.h>
#include <util/gjp.h>
#include <util/lattice.h>
#include <util/link_buffer.h>
#include <util/list.h>
#include <util/qcdio.h>
#include <omp.h>

#include <vector>
#include <map>
#include <stdio.h>
#include <cassert>
CPS_START_NAMESPACE

const int MATRIX_SIZE = 18;

LinkBuffer::LinkBuffer(Lattice &lat_, int buf_sz_not_used):lat(lat_) {
  cname = "LinkBuffer";
  char * fname = "LinkBuffer";
  VRB.Result(cname, fname, "create link buffer object.\n");
  VRB.Func(cname, fname);
  bufferred_gauge_field = NULL;
  bufferred_volume = 0;
  record_bufferred_indexes = false;
}

/*!
  Prints to \c stdout the number of times a lenk has been fetched from the
  buffer and the the number of times that link was on this node.
*/
LinkBuffer::~LinkBuffer(){
  char * fname = "~LinkBuffer()";
  VRB.Result(cname, fname, "destroy link buffer object.\n");
  VRB.Func(cname, fname);
  ClearAll();
}

void Lattice::CreateBuffer(int thickness, const std::vector<int> & indexes) {
  if(LinkBufferIsEnabled()) {
    link_buffer->CreateBuffer(thickness, indexes);
  }
}

void Lattice::CreateBuffer(int thickness, bool record) {
  if(LinkBufferIsEnabled()) {
    link_buffer->CreateBuffer(thickness, record);
  }
}

void Lattice::CreateBuffer(int thickness, LinkBufferShape shape, bool record) {
  if(LinkBufferIsEnabled()) {
    link_buffer->CreateBuffer(thickness, shape, record);
  }
}

void LinkBuffer::CreateBuffer(int thickness, const std::vector<int> & indexes) {
  using namespace std;
  char * fname = "CreateBuffer(indexes)";
  //VRB.Result(cname, fname, "communication size = %d.\n", indexes.size());
  Float dtime = -dclock();
  
  buffer_thickness = thickness;
  record_bufferred_indexes = false;
  bufferred_indexes.clear();
  
  int thick = buffer_thickness;
  int * sites = local_node_sites;
  int * offset = bufferred_dir_offset;
  
  //compute size of buffered volume
  int vol = 1;
  for (int i = 0; i < 4; i++) {
    sites[i] = lat.node_sites[i];
    offset[i] = vol;
    vol *= sites[i] + 2 * thick;
  }
  
  //if the buffered volume has changed, then free the memory
  //for the old buffered volume
  if (vol != bufferred_volume) {
    ClearAll();
  }
  
  //if necessary, allocate the buffered volume
  if (NULL == bufferred_gauge_field) {
    VRB.Result(cname, fname, "malloc buffer vol = %d.\n", vol);
    bufferred_gauge_field = (Matrix *) smalloc(vol * 4 * MATRIX_SIZE * sizeof(Float));
    bufferred_volume = vol;
  }
  
  //the key of send map is a 4-element vector giving the offset of the node to
  //send to. The value is a vector of all the Matrices to send.
  map<vector<int>, vector<Matrix> > sendmap;
  map<vector<int>, int > sendmap_consume;
  vector<int> key(4, 0), pos(4, 0), lpos(4, 0);

  for (int i = 0; i < indexes.size(); i++) {
    GetCoordinateFromIndex(pos.data(), indexes[i]);
    for (int i = 0; i < 4; i++) {
      lpos[i] = pos[i] % sites[i];
      key[i] = pos[i] / sites[i];
      if (lpos[i] < 0) {
        lpos[i] += sites[i];
        key[i]--;
      }
    }
    int mu = indexes[i] % 4;
    sendmap[key].push_back(*lat.GetLink(lpos.data(), mu));
  }
  Float comtime = -dclock();
  map<vector<int>, vector<Matrix> >::iterator it;
  vector<Matrix> recv_vec;
  for (it = sendmap.begin(); it != sendmap.end(); it++) {
    key = it->first;
    //VRB.Result(cname, fname, "communicating : %d %d %d %d\n",
    //    key[0], key[1], key[2], key[3]);
    vector<Matrix> & send_vec = it->second;
    int size = send_vec.size();
    size_t size_bytes = size * MATRIX_SIZE * sizeof(IFloat);
    recv_vec.resize(max(2500, size));
    Matrix * send = send_vec.data();
    Matrix * recv = recv_vec.data();
    for (int i = 0; i < 4; i++) {
      int dis = key[i];
      if (dis < 0) {
        while (dis != 0) {
          getMinusData((IFloat *)recv, (IFloat *)send, size * MATRIX_SIZE, i);
          memcpy(send, recv, size_bytes);
          dis++;
        }
      } else if (dis > 0) {
        while (dis != 0) {
          getPlusData((IFloat *)recv, (IFloat *)send, size * MATRIX_SIZE, i);
          memcpy(send, recv, size_bytes);
          dis--;
        }
      }
    }
    sendmap_consume[key] = 0;
  }
  comtime += dclock();
  //print_flops(cname, fname, 0, comtime);
  for (int i = 0; i < indexes.size(); i++) {
    int index = indexes[i];
    GetCoordinateFromIndex(pos.data(), index);
    for (int i = 0; i < 4; i++) {
      lpos[i] = pos[i] % sites[i];
      key[i] = pos[i] / sites[i];
      if (lpos[i] < 0) {
        lpos[i] += sites[i];
        key[i]--;
      }
    }
    int consume = sendmap_consume[key];
    bufferred_gauge_field[index] = sendmap[key][consume];
    sendmap_consume[key] = consume + 1;

    /* Debug code
    int mu = index % 4;
    Matrix diff_matrix = *lat.GetLink(pos.data(), mu) - bufferred_gauge_field[index];
    if (diff_matrix.ReTr() != 0) {
      VRB.Result(cname, fname, "Mismatch: %d %d %d %d -- %d %d %d %d\n",
          key[0], key[1], key[2], key[3],
          pos[0], pos[1], pos[2], pos[3]);
    }
    // */

  }
  dtime += dclock();
  //print_flops(cname, fname, 0, dtime);
}

void LinkBuffer::CreateBuffer(int thickness, bool record) {
  using namespace std;
  char * fname = "CreateBuffer";
  Float dtime = -dclock();

  buffer_thickness = thickness;
  record_bufferred_indexes = record;
  bufferred_indexes.clear();

  int thick = buffer_thickness;
  int * sites = local_node_sites;
  int * offset = bufferred_dir_offset;

  //compute the size of the link buffer
  int vol = 1;
  for (int i = 0; i < 4; i++) {
    sites[i] = lat.node_sites[i];
    offset[i] = vol;
    vol *= sites[i] + 2 * thick;
  }

  //if the buffered volume has changed, free the old buffer
  if (vol != bufferred_volume) {
    ClearAll();
  }

  //if necessary, allocate memory for a new buffer
  if (NULL == bufferred_gauge_field) {
    VRB.Result(cname, fname, "malloc buffer vol = %d.\n", vol);
    bufferred_gauge_field = (Matrix *) smalloc(vol * 4 * MATRIX_SIZE * sizeof(Float));
    bufferred_volume = vol;
  }

  //the key of send map is a 4-element vector giving the offset of the node to
  //receive from. The value is a vector of all the Matrices to send.
  map<vector<int>, vector<Matrix> > sendmap;
  map<vector<int>, int > sendmap_consume;
  vector<int> pos(4, 0);  //coordinates position of a site relative to this node
  vector<int> lpos(4, 0); //coordinates of a site relative to its home node
  vector<int> key(4, 0);  //offset of a site's home node from this node

  //Populate sendmap with the data that we need to send to other nodes
  for (int x = - thick; x < sites[0] + thick; x++)
    for (int y = - thick; y < sites[1] + thick; y++)
      for (int z = - thick; z < sites[2] + thick; z++)
        for (int t = - thick; t < sites[3] + thick; t++) {
          pos[0] = x;
          pos[1] = y;
          pos[2] = z;
          pos[3] = t;
          for (int i = 0; i < 4; i++) {
            lpos[i] = pos[i] % sites[i];
            key[i] = pos[i] / sites[i];
            if (lpos[i] < 0) {
              lpos[i] += sites[i];
              key[i]--;
            }
          }
          //if (key[0] != 0 || key[1] != 0 || key[2] != 0 || key[3] != 0) {
            vector<Matrix> & vec = sendmap[key];
            for (int mu = 0; mu < 4; mu++) {
              vec.push_back(*lat.GetLink(lpos.data(), mu));
            }
          //}
        }

  Float comtime = -dclock();
  
  vector<Matrix> recv_vec; //will store data received from other nodes

  //Iterate over all the nodes to which we need to send data.
  //We ultimately copy the received data into the corresponding
  //value of sendmap.
  map<vector<int>, vector<Matrix> >::iterator it;
  for (it = sendmap.begin(); it != sendmap.end(); it++) {
    key = it->first; //offset of node
    vector<Matrix> & send_vec = it->second;

    int size = send_vec.size();
    size_t size_bytes = size * MATRIX_SIZE * sizeof(IFloat);
    recv_vec.resize(max(2500, size)); //Ask Luchang why max()

    //Get data from the node with offset given by key. 
    //This may take several rounds of communications.
    Matrix * send = send_vec.data();
    Matrix * recv = recv_vec.data();
    for (int i = 0; i < 4; i++) {
      int dis = key[i];
      if (dis < 0) {
        while (dis != 0) {
          getMinusData((IFloat *)recv, (IFloat *)send, size * MATRIX_SIZE, i);
          memcpy(send, recv, size_bytes);
          dis++;
        }
      } else if (dis > 0) {
        while (dis != 0) {
          getPlusData((IFloat *)recv, (IFloat *)send, size * MATRIX_SIZE, i);
          memcpy(send, recv, size_bytes);
          dis--;
        }
      }
    }
    sendmap_consume[key] = 0;
  }
  comtime += dclock();
  print_flops(cname, fname, 0, comtime);

  //now sendmap[key] is the vector of data recieved from the node
  //pointed to by key.

  //Copy the received data into bufferred_gauge_field
  for (int x = - thick; x < sites[0] + thick; x++)
    for (int y = - thick; y < sites[1] + thick; y++)
      for (int z = - thick; z < sites[2] + thick; z++)
        for (int t = - thick; t < sites[3] + thick; t++) {
          pos[0] = x;
          pos[1] = y;
          pos[2] = z;
          pos[3] = t;
          for (int i = 0; i < 4; i++) {
            lpos[i] = pos[i] % sites[i];
            key[i] = pos[i] / sites[i];
            if (lpos[i] < 0) {
              lpos[i] += sites[i];
              key[i]--;
            }
          }
          int index = BufferredSiteOffset(pos.data()) * 4;
          //if (key[0] != 0 || key[1] != 0 || key[2] != 0 || key[3] != 0) {
            //sendmap_consume[key] keeps track of our progress in consuming the
            //received data in sendmap[key], so that we know which index of 
            //sendmap[key] corresponds to which site.
            int consume = sendmap_consume[key]; 
            vector<Matrix> & vec = sendmap[key];
            for (int mu = 0; mu < 4; mu++) {
              bufferred_gauge_field[index+mu] = vec[consume];
              consume++;

              /* Debug code
              Matrix diff_matrix = *lat.GetLink(pos.data(), mu) - bufferred_gauge_field[index+mu];
              if (diff_matrix.ReTr() != 0) {
                VRB.Result(cname, fname, "Mismatch: %d %d %d %d -- %d %d %d %d\n",
                    key[0], key[1], key[2], key[3],
                    pos[0], pos[1], pos[2], pos[3]);
              }
              // */

            }
            sendmap_consume[key] = consume;
          //}

        }
  dtime += dclock();
  print_flops(cname, fname, 0, dtime);
}


//Computes the offset of a site in a generic subvolume
static inline int GsiteOffset(const int * x,
                              const int * size,
                              const int * origin)
{
  int ret = 4 * (x[0] - origin[0] + size[0] *
                (x[1] - origin[1] + size[1] *
                (x[2] - origin[2] + size[2] *
                (x[3] - origin[3]))));
  return ret;
}

//Computes the coordinates of a site in a generic
//subvolume from the site's index.
static inline void GetCoordinatesFromIndex(int x[4],
                                           const int index,
                                           const int size[4],
                                           const int origin[4])
{
    int j = index;
    x[0] = (j % size[0]) + origin[0]; j /= size[0];
    x[1] = (j % size[1]) + origin[1]; j /= size[1];
    x[2] = (j % size[2]) + origin[2]; j /= size[2];
    x[3] = (j % size[3]) + origin[3];
}

//Copies part of one subvolume into another subvolume. 
//dest_size and dest_origin specify the size and position of the destination subvolume.
//src_size and src_origin specify the size and position of the source subvolume.
//copy_size and copy_origin specify the size and position of the
//subvolume that will actually be copied.
static void CopySubvolume(Matrix* dest,
                          const int dest_size[4],
                          const int dest_origin[4],
                          const Matrix* src,
                          const int src_size[4],
                          const int src_origin[4],
                          const int copy_size[4],
                          const int copy_origin[4])
{
  const int site_bytes = 4 * sizeof(Matrix);
  const int copy_volume = copy_size[0] * copy_size[1] * copy_size[2] * copy_size[3];

  omp_set_num_threads(64);
#pragma omp parallel for
  for(int i = 0; i < copy_volume; i++) {
    int x[4];
    GetCoordinatesFromIndex(x, i, copy_size, copy_origin);

    int dest_offset = GsiteOffset(x, dest_size, dest_origin);
    int src_offset = GsiteOffset(x, src_size, src_origin);

    memcpy(dest + dest_offset, src + src_offset, site_bytes);
  }
}


void LinkBuffer::CreateBuffer(int thickness, LinkBufferShape shape, bool record) {
  char * fname = "CreateBuffer";
  if(shape == BUFFER_SHAPE_HYPERCUBE) {
    CreateBuffer(thickness, record);
    return;
  }

  buffer_thickness = thickness;
  record_bufferred_indexes = record;
  bufferred_indexes.clear();

  int thick = buffer_thickness;
  int * sites = local_node_sites;
  int * offset = bufferred_dir_offset;

  //compute the size of the link buffer
  int vol = 1;
  for (int i = 0; i < 4; i++) {
    sites[i] = lat.node_sites[i];
    offset[i] = vol;
    vol *= sites[i] + 2*thick;
  }

  //if the buffered volume has changed, free the old buffer
  if (vol != bufferred_volume) {
    ClearAll();
  }

  //if necessary, allocate memory for a new buffer
  if (NULL == bufferred_gauge_field) {
    VRB.Result(cname, fname, "malloc buffer vol = %d.\n", vol);
    bufferred_gauge_field = (Matrix *) smalloc(vol * 4 * MATRIX_SIZE * sizeof(Float));
    bufferred_volume = vol;
  }

  const int buffer_size[4] = {sites[0] + 2 * thick,
                              sites[1] + 2 * thick,
                              sites[2] + 2 * thick,
                              sites[3] + 2 * thick};
  const int buffer_origin[4] = {-thick, -thick, -thick, -thick};

  Matrix * local_gauge_field = lat.GaugeField();
  const int local_size[4] = {sites[0], sites[1], sites[2], sites[3]};
  const int local_origin[4] = {0, 0, 0, 0};

  //First do the easy part: copy the local volume into the buffer
  CopySubvolume(bufferred_gauge_field, buffer_size, buffer_origin,
                local_gauge_field, local_size, local_origin,
                local_size, local_origin);

  //Next grab faces from the neighboring nodes
  for(int mu = 0; mu < 4; mu++) {
    assert(local_size[mu] >= thick);

    //Allocate space to store the transmitted slabs
    int slab_size[4] = {local_size[0], local_size[1], local_size[2], local_size[3]};
    slab_size[mu] = thick;
    int slab_vol = slab_size[0] * slab_size[1] * slab_size[2] * slab_size[3];
    
    int slab_bytes = 4 * sizeof(Matrix) * slab_vol;
    int slab_floats = slab_bytes / sizeof(Float);
    Matrix * slab_send = (Matrix *)smalloc(slab_bytes, "slab_send");
    Matrix * slab_recv = (Matrix *)smalloc(slab_bytes, "slab_recv");

    int send_origin[4] = {0, 0, 0, 0};
    int recv_origin[4] = {0, 0, 0, 0};

    for(int send_forward = 0; send_forward <= 1; send_forward++) {
      if(send_forward) { //we will transmit data forward
        send_origin[mu] = local_size[mu] - thick;
        recv_origin[mu] = -thick;
      } else { //we will transmit data backward
        send_origin[mu] = 0;
        recv_origin[mu] = local_size[mu];
      }

      //Assemble the slab that will be transmitted
      CopySubvolume(slab_send, slab_size, send_origin,
                    local_gauge_field, local_size, local_origin,
                    slab_size, send_origin);

      //Transmit the slab
      if(send_forward) getMinusData((IFloat*)slab_recv, (IFloat*)slab_send, slab_floats, mu);
      else              getPlusData((IFloat*)slab_recv, (IFloat*)slab_send, slab_floats, mu);

      //Copy the received slab into the right place in the buffer
      CopySubvolume(bufferred_gauge_field, buffer_size, buffer_origin,
                    slab_recv, slab_size, recv_origin,
                    slab_size, recv_origin);
    }

    sfree(slab_send, "slab_send");
    sfree(slab_recv, "slab_recv");
  }

  //If all we had to buffer was the faces then we are done
  if(shape == BUFFER_SHAPE_LOCAL_AND_FACES) return;

  //Finally we buffer the corners that are two hops from this node
  //(but not the corners that are three or four hops from this node).
  assert(shape == BUFFER_SHAPE_LOCAL_AND_FACES_AND_1ST_CORNERS);
    
  for(int mu = 0; mu < 4; mu++) {
    for(int nu = mu+1; nu < 4; nu++) {
      //Allocate space to store the transmitted corners
      int corner_size[4] = {local_size[0], local_size[1], local_size[2], local_size[3]};
      corner_size[mu] = thick;
      corner_size[nu] = thick;
      int corner_vol = corner_size[0] * corner_size[1] * corner_size[2] * corner_size[3];

      int corner_bytes = 4 * sizeof(Matrix) * corner_vol;
      int corner_floats = corner_bytes / sizeof(Float);
      Matrix * corner_send = (Matrix *)smalloc(corner_bytes, "corner_send");
      Matrix * corner_itmd = (Matrix *)smalloc(corner_bytes, "corner_itmd");
      Matrix * corner_recv = (Matrix *)smalloc(corner_bytes, "corner_recv");

      int send_origin[4] = {0, 0, 0, 0};
      int recv_origin[4] = {0, 0, 0, 0};

      for(int mu_forward = 0; mu_forward <= 1; mu_forward++) {
        for(int nu_forward = 0; nu_forward <= 1; nu_forward++) {
          if(mu_forward) { 
            send_origin[mu] = local_size[mu] - thick;
            recv_origin[mu] = -thick;
          } else { 
            send_origin[mu] = 0;
            recv_origin[mu] = local_size[mu];
          }
          if(nu_forward) { 
            send_origin[nu] = local_size[nu] - thick;
            recv_origin[nu] = -thick;
          } else { 
            send_origin[nu] = 0;
            recv_origin[nu] = local_size[nu];
          }

          //Assemble the corner that will be transmitted
          CopySubvolume(corner_send, corner_size, send_origin,
                        local_gauge_field, local_size, local_origin,
                        corner_size, send_origin);

          //Transmit the corner; this time this requires two hops
          if(mu_forward) getMinusData((IFloat*)corner_itmd, (IFloat*)corner_send, corner_floats, mu);
          else            getPlusData((IFloat*)corner_itmd, (IFloat*)corner_send, corner_floats, mu);

          if(nu_forward) getMinusData((IFloat*)corner_recv, (IFloat*)corner_itmd, corner_floats, nu);
          else            getPlusData((IFloat*)corner_recv, (IFloat*)corner_itmd, corner_floats, nu);

          //Copy the received corner into the right place in the buffer
          CopySubvolume(bufferred_gauge_field, buffer_size, buffer_origin,
                        corner_recv, corner_size, recv_origin,
                        corner_size, recv_origin);
        }
      }

      sfree(corner_send, "corner_send");
      sfree(corner_itmd, "corner_itmd");
      sfree(corner_recv, "corner_recv");
    }
  }
}


std::vector<int> Lattice::GetBufferredIndexes() {
  char * fname = "GetBufferredIndexes";
  if (NULL == link_buffer) {
    ERR.General(cname, fname, "buffer not enabled.\n");
  } else {
    return link_buffer->GetBufferredIndexes();
  }
}

std::vector<int> LinkBuffer::GetBufferredIndexes() {
  using namespace std;
  int size = bufferred_indexes.size();
  vector<int> indexes;
  indexes.resize(size);
  set<int>::iterator it;
  int k = 0;
  for (it = bufferred_indexes.begin(); it != bufferred_indexes.end(); it++) {
    indexes[k] = *it;
    k++;
  }
  return indexes;
}

//----------------------------------------------------------------------
/*!
  Looks for the link \a U_mu(x) in the buffer. If it is not there it is
  brought in.
  \param x The lattice coordinates of the link.
  \param mu The direction index of the link.
  \return The link \a U_mu(x).
*/
//----------------------------------------------------------------------
const Matrix * LinkBuffer::GetBufferedLink(const int *x, int mu) {
  char * fname = "GetBufferedLink()";
  VRB.Func(cname, fname);
  int index = BufferredSiteOffset(x) * 4 + mu;

  //add check to see whether x is in the buffered volume?

  if (index < bufferred_volume * 4) {
    if (record_bufferred_indexes) {
      bufferred_indexes.insert(index);
    }
    return bufferred_gauge_field + index;
  } else {
    //remove this later:
    ERR.General(cname, fname, "Link not in buffer.\n");
    return lat.GetLink(x, mu);
  }
}

//----------------------------------------------------------------------
/*!
  Remove from the buffer all the links with local lattice site \a x
  and direction \a mu.
  \param x The lattice coordinates.
  \param mu The direction index.
*/
//----------------------------------------------------------------------

void LinkBuffer::ClearAll(){
  char * fname = "ClearAll()";
  VRB.Func(cname, fname);
  bufferred_indexes.clear();
  if (NULL != bufferred_gauge_field) {
    VRB.Result(cname, fname, "sfree buffer vol = %d.\n", bufferred_volume);
    sfree(bufferred_gauge_field);
    bufferred_gauge_field = NULL;
    bufferred_volume = 0;
  }
}

void Lattice::ClearAllBufferedLink(){
  if(LinkBufferIsEnabled())
    link_buffer->ClearAll();
}

/*!
  This function does not create a buffer if there already is one.
  \param buf_sz The size of the link buffer.
  \return True if there is a link buffer, false otherwise.
*/
int Lattice::EnableLinkBuffer(int buf_sz){
  //char * fname = "EnableLinkBuffer()";
  //VRB.Func(cname, fname);
  if(link_buffer==0)
    link_buffer = new LinkBuffer(*this, buf_sz);
  if(link_buffer) return 1;
  else return 0;
}

void Lattice::DisableLinkBuffer(){
  if(link_buffer != 0){
    delete link_buffer;
    link_buffer = 0;
  }
}

/*!
  Looks for the link \a U_mu(x) in the buffer. If it is not there it is
  brought in.
  \param x The lattice coordinates of the link.
  \param mu The direction index of the link.
  \return The link \a U_mu(x).
*/
const Matrix * Lattice::GetBufferedLink(const int *x , int mu) {
  //if (IsOnNode(x)) {
  //  /* Debug code
  //  VRB.Result(cname, "", "Read %d %d %d %d mu %d\n", x[0], x[1], x[2], x[3], mu);
  //  Matrix diff_matrix = *(gauge_field+GsiteOffset((int*)x) + mu) - *(link_buffer->GetBufferedLink(x, mu));
  //  if (diff_matrix.ReTr() != 0) {
  //    VRB.Result(cname, "", "Mismatch %d %d %d %d mu %d\n", x[0], x[1], x[2], x[3], mu);
  //  }
  //  // */
  //  return gauge_field+GsiteOffset((int*)x) + mu;
  //}
  if(LinkBufferIsEnabled()) return link_buffer->GetBufferedLink(x, mu);
  //remove this later:
  ERR.General("Lattice", "GetBufferedLink", "Link buffer not enabled\n");
  return GetLink(x,mu);
}

/*!
  \param x The lattice site coordinates
  \return 1 if x is on this node, 0 otherwise.
 */
int Lattice::
IsOnNode(const int * x){
  for(int i=0;i<4;i++)
   if(x[i]<0 || x[i]>=node_sites[i]) return 0; 
  return 1;
}


//------------------------------------------------------------------
//Lattice::BufferedStaple()
/*!
  The staple sum around the link \f$ U_\mu(x) \f$is
  
\f[
  \sum_{v \neq \mu} [                                  
              U_\nu(x+\mu) U_\mu(x+\nu) U^\dagger_\nu(x)                  
           +  U^\dagger_\nu(x+\mu-\nu) U^\dagger_u(x-\nu) U_\nu(x-\nu) ]
\f]

  \param x The coordinates of the lattice site 
  \param u The link direction
  \param stap The computed staple sum.
*/
//------------------------------------------------------------------
void Lattice::
BufferedStaple(Matrix& stap, const int *x, int u){
  //char * fname = "BufferedStaple()";
  //VRB.Func(cname, fname);
  int link_site[4];

  const unsigned CBUF_MODE4 = 0xcca52112;
  const unsigned CBUF_MODE2 = 0xcca52112;

  setCbufCntrlReg(2, CBUF_MODE2);
  setCbufCntrlReg(4, CBUF_MODE4);

  Matrix accumulate;
  accumulate.ZeroMatrix();
  
  for(int i=0;i<4;i++)link_site[i]=x[i]; 

  link_site[u]++;
  
  int dir[3];
  stap.ZeroMatrix();
  
  int u1= u+4;
  for(int v=0; v<8; v++){
    if((v&3)==u) continue;
    dir[0]=v;     //v0
    dir[1]=u1;
    dir[2]=(v+4)&7; //v1
    PathOrdProdPlus(accumulate, link_site, dir, 3);
    //vecAddEquVec((IFloat*) &stap, (IFloat*) &m_tmp,
    //		   MATRIX_SIZE*sizeof(IFloat));
  } 
  moveMem((IFloat*)&stap,  (IFloat*) &accumulate, MATRIX_SIZE*sizeof(IFloat)); 
}


//------------------------------------------------------------------
//Lattice::BufferedRectStaple()
/*! The 5-link rectangle staple sum around the link U_\mu(x) is:

\f[
 \sum_{\nu \neq \mu} \left[\right.                                     
 U_\mu(x+\mu)  U_\nu(x+2\mu)  U^\dagger_\mu(x+\mu+\nu)
 U^\dagger_\mu(x+\nu)  U^\dagger_\nu(x)       \f]\f[
 + U_\mu(x+\mu)  U^\dagger_\nu(x+2\mu-\nu) U^\dagger_\mu(x+\mu-\nu)
 U^\dagger_\mu(x-\nu)  U_\nu(x-\nu)  \f]\f[
 + U_\nu(x+\mu)  U^\dagger_\mu(x+\nu)  U^\dagger_\mu(x-\mu+\nu)
 U^\dagger_\nu(x-\mu) U_\mu(x-\mu)      \f]\f[
 + U^\dagger_\nu(x+\mu-\nu) U^\dagger_\mu(x-\nu)  U^\dagger_\mu(x-\mu-\nu)
 U_\nu(x-\mu-\nu) U_\mu(x-\mu)  \f]\f[
 + U_\nu(x+\mu)  U_\nu(x+\mu+\nu)  U^\dagger_\mu(x+2\nu)
 U^\dagger_\nu(x+\nu)  U^\dagger_\nu(x) \f]\f[
 + U^\dagger_\nu(x+\mu-\nu) U^\dagger_\nu(x+\mu-2\nu) U^\dagger_\mu(x-2\nu)
 U_\nu(x-2\nu)  U_\nu(x-\nu)
 \left.\right]
\f]

  \param x The coordinates of the lattice site 
  \param u The link direction
  \param stap The computed staple sum.
*/
//------------------------------------------------------------------
void Lattice::
BufferedRectStaple(Matrix& stap, const int *x, int u){
  //char * fname = "BufferedRectStaple()";
  //VRB.Func(cname, fname);
  
  int link_site[4];
  const unsigned CBUF_MODE4 = 0xcca52112;
  const unsigned CBUF_MODE2 = 0xcca52112;

  setCbufCntrlReg(2, CBUF_MODE2);
  setCbufCntrlReg(4, CBUF_MODE4);

  Matrix accumulate;
  accumulate.ZeroMatrix();
  
  for(int i=0;i<4;i++)link_site[i]=x[i]; 

  link_site[u]++;
  
  int dir[5];
  stap.ZeroMatrix();
  int u1= u+4;
  for(int v=0; v<8; v++){
    if((v&3)==u) continue;
      
    int v1= (v+4)&7;
      
    dir[0] = v; 
    dir[1] = v;
    dir[2] = u1;
    dir[3] = v1;
    dir[4] = v1;
      
    PathOrdProdPlus(accumulate, link_site, dir, 5); 
    //vecAddEquVec((IFloat*) &stap, (IFloat*) &m_tmp, 
    //		   MATRIX_SIZE*sizeof(IFloat)); 
      
    //dir[0] = v 
    dir[1] = u1;
    //dir[2] = u1;
    //dir[3] = v1; 
    dir[4] = u;
      
    PathOrdProdPlus(accumulate, link_site, dir, 5); 
    //vecAddEquVec((IFloat*) &stap, (IFloat*) &m_tmp, 
    //		   MATRIX_SIZE*sizeof(IFloat)); 
      
    dir[0] = u; 
    dir[1] = v;
    //dir[2] = u1;
    dir[3] = u1; 
    dir[4] = v1;
    PathOrdProdPlus(accumulate, link_site, dir, 5); 
    //vecAddEquVec((IFloat*) &stap, (IFloat*) &m_tmp, 
    //		   MATRIX_SIZE*sizeof(IFloat)); 
    
  }
  moveMem((IFloat*)&stap, (IFloat*)&accumulate, MATRIX_SIZE*sizeof(IFloat));
}


//------------------------------------------------------------------
/*! The chair shaped 5-link staple sum around the link U_\mu(x) is:
  
\f[
 \sum_{\pm \nu, |\nu|\neq \mu} \sum_{\pm \rho, |\rho|\neq \nu, |\rho|\neq \mu}(
\left[\right.
 U_\rho(x+\mu) U_\nu(x+\mu+\rho) U_{-\rho}(x+\mu+\nu+\rho)
 U_{-\mu}(x+\mu+\nu) U_{-\nu}(x+\nu)   \f]\f[
 + U_\rho(x+\mu) U_\nu(x+\mu+\rho) U_{-\mu}(x+\mu+\nu+\rho)
 U_{-\nu}(x+\nu+\rho) U_{-\rho}(x+\rho)    \f]\f[
 + U_\rho(x+\mu) U_{-\mu}(x+\mu+\rho) U_\nu(x+\rho)
 U_{-\rho}(x+\rho+\nu) U_{-\nu}(x+\nu)        
 \left.\right]
\f]

  \param x The coordinates of the lattice site 
  \param u The link direction
  \param stap The computed staple sum.
*/
//------------------------------------------------------------------
void Lattice::
BufferedChairStaple(Matrix& stap, const int *x, int u){
 
  int link_site[4];
  
  for(int i=0;i<4;i++) link_site[i]=x[i]; 

  link_site[u]++;

  const unsigned CBUF_MODE4 = 0xcca52112;
  const unsigned CBUF_MODE2 = 0xcca52112;

  setCbufCntrlReg(2, CBUF_MODE2);
  setCbufCntrlReg(4, CBUF_MODE4);

  Matrix accumulate;
  accumulate.ZeroMatrix();
  
  int dir[5];
  Matrix m_tmp;
  
  int u1= u+4;
  
  for(int v=0; v<8; v++){
    
    if((v&3)==u) continue;
      
    int v1= (v+4)&7;
      
    for(int w=0; w<8; w++){
      if((w&3) ==u || (w&3) == (v&3) ) continue;
	  
      int w1= (w+4)&7;
	  
      dir[0] = w; 
      dir[1] = v;
      dir[2] = w1;
      dir[3] = u1;
      dir[4] = v1;

      PathOrdProdPlus(accumulate, link_site, dir, 5); 
      //vecAddEquVec((IFloat*) &stap, (IFloat*) &m_tmp, 
      //          MATRIX_SIZE*sizeof(IFloat)); 
	  
      //dir[0] = w 
      //dir[1] = v
      dir[2] = u1;
      dir[3] = v1; 
      dir[4] = w1;
      PathOrdProdPlus(accumulate, link_site, dir, 5); 
      //vecAddEquVec((IFloat*) &stap, (IFloat*) &m_tmp, 
      //            MATRIX_SIZE*sizeof(IFloat)); 
	  
      //dir[0] = w; 
      dir[1] = u1;
      dir[2] = v;
      dir[3] = w1; 
      dir[4] = v1;
      PathOrdProdPlus(accumulate, link_site, dir, 5); 
      //vecAddEquVec((IFloat*) &stap, (IFloat*) &m_tmp, 
      //             MATRIX_SIZE*sizeof(IFloat)); 

    }
  }
  moveMem((IFloat *) &stap, (IFloat*)&accumulate, MATRIX_SIZE*sizeof(IFloat));
}
//------------------------------------------------------------------------
/*!
  The staple sum around the link U_\mu(x) is
\f[
\sum_{\pm \nu, |\nu|\neq \mu} \sum_{\pm \rho, |\rho|\neq \nu, |\rho|\neq \mu}

     U_\nu(x+\mu) U_\rho(x+\mu+\nu) U_{-\mu}(x+\mu+\nu+\rho) U_{-\nu}(x+\nu+\rho) U_{-\rho}(x+\rho)
    
\f]

  \param x The coordinates of the lattice site 
  \param u The link direction
  \param stap The computed staple sum.
*/     
//------------------------------------------------------------------------
void Lattice:: 
BufferedCubeStaple(Matrix &stap, const int *x, int u){ 
  int link_site[4];

  const unsigned CBUF_MODE4 = 0xcca52112;
  const unsigned CBUF_MODE2 = 0xcca52112;

  setCbufCntrlReg(2, CBUF_MODE2);
  setCbufCntrlReg(4, CBUF_MODE4);

  Matrix accumulate;
  accumulate.ZeroMatrix();

  for(int i=0;i<4;i++) link_site[i]=x[i]; 

  link_site[u]++;
  
  int dir[5];
  stap.ZeroMatrix();
  Matrix m_tmp;
  for(int v=0; v<8; v++){
    if((v&3)==u) continue;
    int v1 = (v+4)&7;
    for(int w=0; w<8; w++){
      if((w&3) ==u || (w&3) == (v&3) ) continue;
      dir[0] = v; 
      dir[1] = w;
      dir[2] = u+4;
      dir[3] = v1;
      dir[4] = (w+4)&7;

      PathOrdProdPlus(accumulate, link_site, dir, 5); 
      //vecAddEquVec((IFloat*) &stap, (IFloat*) &m_tmp, 
      //	       MATRIX_SIZE*sizeof(IFloat)); 
      
    }
  }
  moveMem((IFloat*)&stap, (IFloat*)&accumulate, MATRIX_SIZE*sizeof(IFloat));
}


//-------------------------------------------------------------------
/*!
  Given the starting site x, the directions of each step on the path
  and the number of steps. calculate the path ordered product of
  all the links along the path and add it to the given matrix m.
  Each direction is one of 0, 1, 2, 3, 4, 5, 6 or 7} corresponding to
  the directions X, Y, Z, T, -X, -Y, -Z and -T respectively.

  \param m The initial matrix.
  \param x The coordinates of the starting point of the path
  \param dirs The list of directions.
  \param n The number of links in the path.
  \post \a The product along the path is added to \a m.
*/
// 
// in this implementation, each link is retrieved from other sites,
// and assembled on the local node.
//
// in another implementation, one could imagine passing a matrix
// around along the path, and multiply the links to it along the way.
// this would save a lot of communication, in
// many cases, but not optimal if there is only one offnode link.
// and the link buffer would not be of much use in this case.
//-------------------------------------------------------------------
void Lattice::
PathOrdProdPlus(Matrix & mat, const int * x, const int* dirs, int n){
  //char * fname = "PathOrdProd"; 
  //VRB.Flow(cname, fname,"(,,,%d)\n",n);

  Matrix result, result1, mat1;
  Matrix * result_mp = &result;
  Matrix * result1_mp = &result1;

  int abs_dir;
  int dir_sign;
  int link_site[4];
  
  const Matrix * p1;
  //Matrix m1, m2, m3;

  int i;
  for(i=0;i<4;i++)link_site[i]=x[i]; 

  //deal with the first link
  //------------------------------
  abs_dir = dirs[0]&3; 
    //abs_dir is {0,1,2,3}
  dir_sign = dirs[0]>>2; 
    //dir_sign is {0, 1}

  link_site[abs_dir] -= dir_sign; 
    //if dir_sign == 1, the link is at x-n_v

  p1 = GetBufferedLink(link_site, abs_dir);  
    //get the first link 

  link_site[abs_dir] += 1-dir_sign; 
    //if dir_sign == 0, march on to the next site, if dir_sign == 1, we have
    //already moved.

  if (dir_sign) 
      result1_mp->Dagger((IFloat*)p1); 
        //if dir_sign==1 the link is going backward so get its dagger
  else 
      moveMem((IFloat*) result1_mp, (IFloat*)p1, 
              MATRIX_SIZE * sizeof(IFloat));
        //simply move to the cram
  
  for(i=1;i<n;i++){
    abs_dir = dirs[i]&3; 
    dir_sign = dirs[i]>>2; 
    
    link_site[abs_dir] -= dir_sign; 

    p1 = GetBufferedLink(link_site, abs_dir);  

    link_site[abs_dir] += 1-dir_sign; 

    //put the next link on the path in mat1
    //--------------------------------------
    if(dir_sign){  
      mat1.Dagger((IFloat*)p1); 
    }
    else 
      moveMem((IFloat*)&mat1, (IFloat*)p1, 
              MATRIX_SIZE*sizeof(IFloat));

    if(i!=n-1)
      mDotMEqual((IFloat*) result_mp, (IFloat*)result1_mp, (IFloat*)&mat1);
        //if not the last link on the path, just multiply to the earlier result
    else
      mDotMPlus((IFloat*) &mat, (IFloat*) result1_mp, (IFloat*) &mat1);
        //if the last link, multiply and add to mat.

    Matrix * tmp_p = result1_mp; result1_mp=result_mp; result_mp = tmp_p;
    //swap result_mp and result1_mp;
  }
}

//-------------------------------------------------------------------
/*!
  Given the starting site x, the directions of each step on the path
  and the number of steps. calculate the path ordered product of
  all the links along the path and add it to the given matrix m.
  Each direction is one of 0, 1, 2, 3, 4, 5, 6 or 7} corresponding to
  the directions X, Y, Z, T, -X, -Y, -Z and -T respectively.

  The idea is whenever the path hits a boundary the current partial result is 
  passed to the next processer on the path, which will calculate the part of 
  the product that is on this node and pass on, so the result ends up on the 
  last processor where the path stops.
  
  \param m The product along the path.
  \param x The coordinates of the starting point of the path
  \param dirs The list of directions.
  \param n The number of links in the path.
*/
//-------------------------------------------------------------------

void Lattice::
PathOrdProd(Matrix & mat, const int * x, const int* dirs, int n){
  //char * fname = "PathOrdProd"; 
  //VRB.Flow(cname, fname,"(,,,%d)\n",n);

  int abs_dir;
  int dir_sign;
  int link_site[4];
 
  const Matrix * p1;
  Matrix m1, m2,  m4;
  Matrix mat1, mat2, mat3;
  Matrix * r1_mp = &mat1;
  Matrix * r_mp  = &mat2;
  Matrix * buf1_mp=&mat3;

  int i;
  for(i=0;i<4;i++) { 
    int l_x = x[i];
    while(l_x<0) l_x += node_sites[i]; 
    while(l_x>=node_sites[i]) l_x -= node_sites[i];
    link_site[i]=l_x;
  }

  //deal with the first link
  //-------------------------
  abs_dir = dirs[0]&3; 
  dir_sign = dirs[0]>>2; 

  link_site[abs_dir] -= dir_sign; 
    //go to the site where the link is stored
  if(link_site[abs_dir]<0) link_site[abs_dir] +=node_sites[abs_dir];
    //if out of boundary move to the other boundary in this direction
    //where a path is entering the territory of this node.

  p1 = gauge_field+GsiteOffset((int*)link_site) + abs_dir;
    //get the fisrt link

  link_site[abs_dir] += 1-dir_sign;
    //march on to the next site, if haven't done so 

  if(link_site[abs_dir]==node_sites[abs_dir]){ 
      //just hit another boundary we need to pass the partial result
      //on to the positive direction and get a link from the minus direction
    link_site[abs_dir]=0;
    getMinusData((IFloat*) &m1, (IFloat*) p1, MATRIX_SIZE, abs_dir); 
    moveMem((IFloat*)buf1_mp, (IFloat*)&m1, MATRIX_SIZE);
      //move it to the cram
  }
  else
    moveMem((IFloat*)buf1_mp, (IFloat*)p1, MATRIX_SIZE);
      //didn't hit the boundary, directly move the cram

#define SWAP(a, b) {Matrix* tmp_mp = (a); (a)=(b); (b)=tmp_mp;}

  if (dir_sign) 
     //link going backward so take its hermite conjugate 
    r1_mp->Dagger((IFloat*) buf1_mp); 
  else 
     //put the most recent result in r1_mp
    SWAP(r1_mp, buf1_mp) 
      

  for(i=1;i<n;i++){
    abs_dir = dirs[i]&3; 
      //abs_dir is {0,1,2,3}
    dir_sign = dirs[i]>>2; 
      //dir_sign is {0, 1}
    
    link_site[abs_dir] -= dir_sign; 
      //go to the site of the link

    if(link_site[abs_dir]<0){
        //the next link is offnode pass the current result to that node
        //and receive one from the plus direction
      link_site[abs_dir] += node_sites[abs_dir];
      moveMem((IFloat*)&m1, (IFloat*)r1_mp, MATRIX_SIZE); 
        //move from cram to dram, so scu can access.
      getPlusData((IFloat*)&m2, (IFloat*)&m1, MATRIX_SIZE, abs_dir); 
      moveMem((IFloat*)r1_mp, (IFloat*)&m2, MATRIX_SIZE); 
        //move from dram to cram to speed matrix multiplication.
    } 
    
    p1 = gauge_field+GsiteOffset((int*)link_site) + abs_dir;
      //get the next link
    
    if(dir_sign){ 
       //going backward, take the hermite conjugate
      buf1_mp->Dagger((IFloat*)p1); 
    }
    else{ 
      moveMem((IFloat*)buf1_mp,(Matrix *)p1, MATRIX_SIZE*sizeof(IFloat));
       //simply move to the cram
    }

    mDotMEqual((IFloat*) r_mp, (IFloat*)r1_mp, (IFloat*)buf1_mp);
      //every thing is in cram, just do the multiplication

    link_site[abs_dir] += 1-dir_sign;
       //march on to the next site 

    if(link_site[abs_dir]==node_sites[abs_dir]){
       //just hit the wall in the plus direction, pass on
      link_site[abs_dir]=0;
      moveMem((IFloat*)&m1, (IFloat*)r_mp, MATRIX_SIZE); 
       //cram to dram
      getMinusData((IFloat*) &m2, (IFloat*) &m1, MATRIX_SIZE, abs_dir);
      moveMem((IFloat*)r1_mp, (IFloat*)&m2, MATRIX_SIZE*sizeof(IFloat)); 
       //dram to cram
    }
    else 
       //didn't hit the wall, swap r1_mp, r_mp so r1_mp will hold the most
       //recent result of the products.
      SWAP(r1_mp, r_mp)
        
#undef SWAP
  }
  moveMem((IFloat*)&mat, (IFloat*)r1_mp, MATRIX_SIZE* sizeof(IFloat));
}


CPS_END_NAMESPACE
