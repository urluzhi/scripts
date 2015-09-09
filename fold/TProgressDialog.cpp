/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Authors: Chris Connett and Andrew Yohn, 2006
 */

#include "TProgressDialog.h"


using namespace std;



const char TProgressDialog::spinchars[] = {'/', '-', '\\', '|'};

TProgressDialog::TProgressDialog(ostream &_s)
  : spinstate(0),
    s(_s) {
}

TProgressDialog::~TProgressDialog() {
}

void TProgressDialog::update(int percent) {

  s << "\r";
  s.width(3);
  s << percent << "% [";

  int i = 0;
  for (int i = 0; i < 100; i += 2) {
    if (i <= percent) {
      s << "=";
    } else {
      s << " ";
    }
  }

  s << "] ";
  if (percent < 100) {
    s << spinchars[spinstate] << "                     ";
    // Exactly 21 spaces puts cursor at last character of an
    // 80-character wide terminal
  } else {
    s << " \n";
  }
  s.flush();
  ++spinstate %= 4;

}
