#include "main_window.h"
#include "ui_main_window.h"
//#include "global.h"

MainWindow::MainWindow(QWidget *parent) :
  QMainWindow(parent),
  ui(new Ui::MainWindow)
{
  ui->setupUi(this);
  auto baseSize = ui->opengl_widget->baseSize();
  ui->opengl_widget->resize(baseSize);
}

MainWindow::~MainWindow()
{
  delete ui;
}

