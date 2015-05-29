#ifndef HAMMOCKSUPPORT_H
#define HAMMOCKSUPPORT_H

template <class Real> class ReflectiveObjectRenderer;
class HammockSupport
{
public:
  HammockSupport();
  void Render();
  ~HammockSupport();
private:
  ReflectiveObjectRenderer<double>* support_;
};

#endif // HAMMOCKSUPPORT_H
