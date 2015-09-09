/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Authors: Chris Connett and Andrew Yohn, 2006
 */

#ifndef TEXTPROGRESSBAR_H
#define TEXTPROGRESSBAR_H

#include <iostream>

class TProgressDialog {
private:
  static const char spinchars[];
  int spinstate;
  std::ostream &s;
public:
  TProgressDialog(std::ostream &_s = std::cout);
  virtual ~TProgressDialog();
  void update(int percent);
};

#endif
