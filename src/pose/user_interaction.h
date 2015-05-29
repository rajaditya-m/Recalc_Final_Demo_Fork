#ifndef USER_INTERACTION_H
#define USER_INTERACTION_H
class CudaSkeleton;
class QMouseEvent;
namespace body {
extern double moving_target[3];
bool HandleMousePress(QMouseEvent * event);
bool HandleMouseMove(QMouseEvent * event);
bool HandleMouseRelease(QMouseEvent * event);
}

namespace cloth {
bool HandleMousePress(QMouseEvent * event);
bool HandleMouseMove(QMouseEvent * event);
bool HandleMouseRelease(QMouseEvent * event);
}
#endif // USER_INTERACTION_H
