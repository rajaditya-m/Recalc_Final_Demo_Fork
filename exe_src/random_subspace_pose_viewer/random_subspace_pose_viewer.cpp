#include <fstream>
#include <QApplication>
#include "main_window.h"
#include "rainbow_color.h"
#include "open_gl_qt.h"
#include "global.h"
#include "string_formatter.h"
#include "multi_domain_tet.h"
#include "input_handler.h"


namespace random_subspace_pose_viewer {
MultiDomainTet* multi_tet = NULL;
int pose_num = 100;
int current_pose = 0;
int render_mode = 0;
std::vector<MultiDomainTet::Vec> pose_q;

void Init() {
  multi_tet = new MultiDomainTet(GetMeshFile());
  multi_tet->LoadPartitionInfo(GetPartitionFolder());
  auto GetSubspaceFileName = [](int part_id) -> const char* {
    static std::string name;
//    name = dj::Format("%s/modal_basis/partition_%d.basis.bin", GetDataFolder(), part_id);
    name = dj::Format("%s/modal_basis/partition_%d.rigid_basis.bin", GetDataFolder(), part_id);
    return name.c_str();
  };
  multi_tet->LoadSubspace(GetSubspaceFileName, MultiDomainTet::kBinary);
  multi_tet->AssembleGlobalBasis();
  // read randomly sampled q
  {
    std::ifstream in(dj::Format("%z/random_pose.txt", GetDataFolder()), std::ios::binary);
    int part_num;
    in >> pose_num >> part_num;
    ASSERT(part_num == multi_tet->part_num_);
    for (int p = 0; p < part_num; ++p) {
      int basis_size;
      in >> basis_size;
      ASSERT(basis_size == multi_tet->part_basis_size_[p]);
    }
    pose_q.resize(pose_num);
    for (int i = 0; i < pose_num; ++i) {
      pose_q[i] = MultiDomainTet::Vec::Zero(multi_tet->total_basis_num_);
      for (int n = 0; n < multi_tet->total_basis_num_; ++n) {
        in >> pose_q[i][n];
      }
    }
  }
  current_pose = 0;
  multi_tet->q_ = pose_q[0];
  multi_tet->UpdatePosition();
}

void NextPose() {
  current_pose = (current_pose + 1) % pose_num;
  multi_tet->q_ = pose_q[current_pose];
  P(multi_tet->q_.norm());
  multi_tet->UpdatePosition();
}
void PrevPose() {
  current_pose--;
  if (current_pose < 0) current_pose = pose_num - 1;
  multi_tet->q_ = pose_q[current_pose];
  P(multi_tet->q_.norm());
  multi_tet->UpdatePosition();
}
}

using namespace random_subspace_pose_viewer;
class RandomSubspacePoseViewer : public InputHandler {
public:
  // Mouse handlers
  virtual int HandleMouseMove(QMouseEvent* e) { Q_UNUSED(e); return kNotHandled; }
  virtual int HandleMousePress(QMouseEvent* e) { Q_UNUSED(e); return kNotHandled; }
  virtual int HandleMouseRelease(QMouseEvent* e) { Q_UNUSED(e); return kNotHandled; }
  // Keyboard handlers
  virtual int HandleKeyPress(QKeyEvent* e) {
    switch (e->key()) {
      case Qt::Key_R:
        render_mode = (render_mode + 1) / MultiDomainTet::kRenderModeNum;
        break;
      case Qt::Key_J:
        PrevPose();
        break;
      case Qt::Key_K:
        NextPose();
        break;
      default:
        break;
    }
    return kNotHandled;
  }

  virtual int HandleKeyRelease(QKeyEvent* e) { Q_UNUSED(e); return kNotHandled; }

  virtual int Render() {
    multi_tet->Render(MultiDomainTet::kRenderMode[render_mode], multi_tet->X);
    QFont sansFont("Helvetica [Cronyx]", 14);
    sansFont.setFamily("sans serif");
    int font_height = QFontMetrics(sansFont).height();
    glColor3fv(kGreen());
    global::gl->renderText(5, font_height + 5,
                           QString("pose id: %1 / %2").arg(current_pose).arg(pose_num - 1),
                           sansFont);
    return kNotHandled;
  }
  virtual int Idle() { return kNotHandled; }
  virtual ~RandomSubspacePoseViewer() {}

private:
  RandomSubspacePoseViewer() {
    Init();
  }
  DECLARE_SINGLETON_CLASS(RandomSubspacePoseViewer);
};

int main(int argc, char *argv[]) {
  QApplication a(argc, argv);
  setlocale(LC_NUMERIC, "C");
  MainWindow w;
    #ifdef WINDOW_TITLE
#if defined(_WIN32) || defined(_WIN64)
  QString title(STRINGIZE_TOKEN(WINDOW_TITLE));
#else
  QString title(WINDOW_TITLE);
#endif
  w.setWindowTitle(title);
#endif
  w.show();
  w.hide();

  global::gl->InstallHandler(RandomSubspacePoseViewer::Instance());
  w.show();
  return a.exec();
}

