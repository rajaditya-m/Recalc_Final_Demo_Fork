#ifndef TRIANGLE_MESH_CREATOR_H
#define TRIANGLE_MESH_CREATOR_H
void CreateCylinder(double radius, double height, int slice, int stack,
                    std::vector<double>& verts, std::vector<double>& normal, std::vector<int>& tri);

void CreateCylinder(float radius, float height, int slice, int stack,
                    std::vector<float>& verts, std::vector<float>& normal, std::vector<int>& tri);

void CreateSphere(double radius, int slice, int stack,
                  std::vector<double>& verts, std::vector<double>& normal, std::vector<int>& tri,
                  std::vector<double>* tex_coord = nullptr);

void CreateSphere(float radius, int slice, int stack,
                  std::vector<float>& verts, std::vector<float>& normal, std::vector<int>& tri,
                  std::vector<float>* tex_coord = nullptr);

//                    std::vector<double>* tex_coord = nullptr);
void ExpandTriangleMesh(std::vector<int>& tris, std::vector<double>& verts, std::vector<double>& normal,
                        std::vector<double>& expanded_verts, std::vector<double>& expanded_normal);
#endif // TRIANGLE_MESH_CREATOR_H
