// vim: set ts=2 sw=2 expandtab:
#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of the LinkBuffer class.

  $Id $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004/08/18 11:57:37 $
//  $Header: /space/cvs/cps/cps++/include/util/link_buffer.h,v 1.4 2004/08/18 11:57:37 zs Exp $
//  $Id: link_buffer.h,v 1.4 2004/08/18 11:57:37 zs Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: link_buffer.h,v $
//  $Revision: 1.4 $
//  $Source: /space/cvs/cps/cps++/include/util/link_buffer.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#ifndef LINK_BUFFER_H
#define LINK_BUFFER_H                    //!< Prevent multiple inclusion.

CPS_END_NAMESPACE
#include <util/list.h>
#include <set>
#include <vector>
CPS_START_NAMESPACE

class Matrix;
class Lattice;

//-----------------------------------------------------------------
//
//! LinkBuffer is a class that buffers the off-node links. 
// 
//-----------------------------------------------------------------

class LinkBuffer{
  char * cname;

  Lattice &lat;
  //reference to the lattice

  int bufferred_volume;
  int buffer_thickness;
  int local_node_sites[4];
  int bufferred_dir_offset[4];

  int BufferredSiteOffset(const int * x) {
    int thick = buffer_thickness;
    int * offsets = bufferred_dir_offset;
    return (x[0] + thick) * offsets[0]
         + (x[1] + thick) * offsets[1]
         + (x[2] + thick) * offsets[2]
         + (x[3] + thick) * offsets[3];
  }

  void GetCoordinateFromIndex(int * x, int index) {
    int thick = buffer_thickness;
    int * offsets = bufferred_dir_offset;
    index /= 4;
    x[3] = index / offsets[3] - thick;
    index %= offsets[3];
    x[2] = index / offsets[2] - thick;
    index %= offsets[2];
    x[1] = index / offsets[1] - thick;
    index %= offsets[1];
    x[0] = index / offsets[0] - thick;
  }

  Matrix * bufferred_gauge_field;

  bool record_bufferred_indexes;
  std::set<int> bufferred_indexes;

 public: 

//! Allocate the link buffer.  
/*!
  Allocates a buffer for \a size links on lattice \a lat
  \param lat The lattice on which the links live.
  \param size : no longer useful, see CreateBuffer.
*/
  LinkBuffer(Lattice & lat, int size); 
	
  ~LinkBuffer();

  void CreateBuffer(int thickness, bool record = false);
  void CreateBuffer(int thickness, const std::vector<int> & indexes);
  void CreateBuffer(int thickness, LinkBufferShape shape, bool record = false);


  std::vector<int> GetBufferredIndexes();

  //----------------------------------------------------------------
  // GetBufferedLink looks into the buffer to search for the a given 
  // link if it's not there, it gets the link from the neighbours 
  //----------------------------------------------------------------
  //! Gets a link from the buffer.
  const Matrix * GetBufferedLink(const int *x, int mu);
  
  //! Removes all links from the buffer.
  void ClearAll();

};

#endif


CPS_END_NAMESPACE
