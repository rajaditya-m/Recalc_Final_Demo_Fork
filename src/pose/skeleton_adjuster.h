#ifndef SKELETON_ADJUSTER_H
#define SKELETON_ADJUSTER_H

#include "input_handler.h"
#include "global.h"

class SkeletonAdjuster : public InputHandler
{
public:
  virtual int HandleKeyPress(QKeyEvent *e)
  {
    bool handled = true;
   const float kTranslateStep = 0.01f;
    switch (e->key()) {
    case Qt::Key_1: {
      skeleton->Translate(-kTranslateStep, 0, 0);
      skeleton->UpdatePosition();
    }
    break;
    case Qt::Key_2: {
      skeleton->Translate(+kTranslateStep, 0, 0);
      skeleton->UpdatePosition();
    }
    break;
    case Qt::Key_3: {
      skeleton->Translate(0, -kTranslateStep, 0);
      skeleton->UpdatePosition();
    }
    break;
    case Qt::Key_4: {
      skeleton->Translate(0, +kTranslateStep, 0);
      skeleton->UpdatePosition();
    }
    break;
    case Qt::Key_5: {
      skeleton->Translate(0, 0, -kTranslateStep);
      skeleton->UpdatePosition();
    }
    break;
    case Qt::Key_6: {
      skeleton->Translate(0, 0, +kTranslateStep);
      skeleton->UpdatePosition();
    }
    break;
    case Qt::Key_7: {
      skeleton->Scale(0.99);
      skeleton->UpdatePosition();
    }
    break;
    case Qt::Key_8: {
      skeleton->Scale(1.01);
      skeleton->UpdatePosition();
    }
    break;
    case Qt::Key_9:
      skeleton->ExportFakeBVH("adjust.fbvh");
      L("Saved");
      break;
    case Qt::Key_T:
      point[0] -= kTranslateStep;
      P(point);
      break;
    case Qt::Key_Y:
      point[0] += kTranslateStep;
      P(point);
      break;
    case Qt::Key_U:
      KK;
      point[1] -= kTranslateStep;
      P(point);
      break;
    case Qt::Key_I:
      point[1] += kTranslateStep;
      P(point);
      break;
    case Qt::Key_O:
      point[2] -= kTranslateStep;
      P(point);
      break;
    case Qt::Key_P:
      point[2] += kTranslateStep;
      P(point);
      break;
    default:
      handled = false;
    }
    return (handled) ? 15 : -1;
  }

private:
  SkeletonAdjuster() : InputHandler(10) {}
  DECLARE_SINGLETON_CLASS(SkeletonAdjuster)
};

#endif // SKELETON_ADJUSTER_H
