#ifndef COIN_SOROTOR_H
#define COIN_SOROTOR_H

/**************************************************************************\
 *
 *  This file is part of the Coin 3D visualization library.
 *  Copyright (C) 1998-2005 by Systems in Motion.  All rights reserved.
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  ("GPL") version 2 as published by the Free Software Foundation.
 *  See the file LICENSE.GPL at the root directory of this source
 *  distribution for additional information about the GNU GPL.
 *
 *  For using Coin with software that can not be combined with the GNU
 *  GPL, and for taking advantage of the additional benefits of our
 *  support services, please contact Systems in Motion about acquiring
 *  a Coin Professional Edition License.
 *
 *  See <URL:http://www.coin3d.org/> for more information.
 *
 *  Systems in Motion, Postboks 1283, Pirsenteret, 7462 Trondheim, NORWAY.
 *  <URL:http://www.sim.no/>.
 *
\**************************************************************************/

#include <Inventor/nodes/SoSubNode.h>
#include <Inventor/nodes/SoRotation.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoSFBool.h>
#include <Inventor/SbVec3f.h>
#include <Inventor/SbTime.h>

class SoSensor;
class SoRotorP;

class COIN_DLL_API SoRotor : public SoRotation {
  typedef SoRotation inherited;

  SO_NODE_HEADER(SoRotor);

public:
  static void initClass(void);
  SoRotor(void);

  SoSFFloat speed;
  SoSFBool on;

protected:
  virtual ~SoRotor();

private:
  SoRotorP * pimpl;
  static void oneshotSensorCB(void * d, SoSensor * s);
  static void fieldSensorCB(void * d, SoSensor * s);
  void setRotation(void);
};

#endif // !COIN_SOROTOR_H
