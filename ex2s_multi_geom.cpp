//                                MFEM Example 2
//
// Compile with: make ex2
//
// Sample runs:  ex2 -m ../data/beam-tri.mesh
//               ex2 -m ../data/beam-quad.mesh
//               ex2 -m ../data/beam-tet.mesh
//               ex2 -m ../data/beam-hex.mesh
//               ex2 -m ../data/beam-wedge.mesh
//               ex2 -m ../data/beam-quad.mesh -o 3 -sc
//               ex2 -m ../data/beam-quad-nurbs.mesh
//               ex2 -m ../data/beam-hex-nurbs.mesh
//
// Description:  This example code solves a simple linear elasticity problem
//               describing a multi-material cantilever beam.
//
//               Specifically, we approximate the weak form of -div(sigma(u))=0
//               where sigma(u)=lambda*div(u)*I+mu*(grad*u+u*grad) is the stress
//               tensor corresponding to displacement field u, and lambda and mu
//               are the material Lame constants. The boundary conditions are
//               u=0 on the fixed part of the boundary with attribute 1, and
//               sigma(u).n=f on the remainder with f being a constant pull down
//               vector on boundary elements with attribute 2, and zero
//               otherwise. The geometry of the domain is assumed to be as
//               follows:
//
//                                 +----------+----------+
//                    boundary --->| material | material |<--- boundary
//                    attribute 1  |    1     |    2     |     attribute 2
//                    (fixed)      +----------+----------+     (pull down)
//
//               The example demonstrates the use of high-order and NURBS vector
//               finite element spaces with the linear elasticity bilinear form,
//               meshes with curved elements, and the definition of piece-wise
//               constant and vector coefficient objects. Static condensation is
//               also illustrated.
//
//               We recommend viewing Example 1 before viewing this example.

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;


void D_exact(const Vector &, Vector &);

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = "../data/beam-tri.mesh";
   int order = 1;
   bool static_cond = false;
   bool visualization = 1;
   const char *geometry = "beam";
   int ref_levels = 0;
   int higher_order = 4;
   const char *file = "xxx_";
   const char *solution_file;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
                  "--no-static-condensation", "Enable static condensation.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&geometry, "-g", "--geometry",
                  "Problem geometry, choose between: beam, plate or seat.");
   args.AddOption(&ref_levels, "-l", "--hlevel",
                  "Level of h refinement.");
   args.AddOption(&file, "-f", "--fileout", "Name of the output file.");
   args.AddOption(&solution_file, "-s2", "--solution2","Grid function for higher order solution.");
   args.AddOption(&higher_order, "-ho", "--higher_order","Finite element order (higher order polynomial degree).");

   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   // 2. Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral or hexahedral elements with the same code.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();

   if (mesh->attributes.Max() < 2 || mesh->bdr_attributes.Max() < 2)
   {
      cerr << "\nInput mesh should have at least two materials and "
           << "two boundary attributes! (See schematic in ex2.cpp)\n"
           << endl;
      return 3;
   }

   // 3. Select the order of the finite element discretization space. For NURBS
   //    meshes, we increase the order by degree elevation.
   if (mesh->NURBSext)
   {
      mesh->DegreeElevate(order, order);
   }

   // 4. Refine the mesh to increase the resolution. In this example we do
   //    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
   //    largest number that gives a final mesh with no more than 5,000
   //    elements.
   {
      cout << "Number of elements: " << mesh->GetNE() << endl;
      cout << "Refinement level: "   << ref_levels    << endl;
      for (int l = 0; l < ref_levels; l++)
      {
         mesh->UniformRefinement();
      }
   }

   // 5. Define a finite element space on the mesh. Here we use vector finite
   //    elements, i.e. dim copies of a scalar finite element space. The vector
   //    dimension is specified by the last argument of the FiniteElementSpace
   //    constructor. For NURBS meshes, we use the (degree elevated) NURBS space
   //    associated with the mesh nodes.
   FiniteElementCollection *fec, *fec_ho;
   FiniteElementSpace *fespace, *fespace_ho;

   if (mesh->NURBSext)
   {
      fec = NULL;
      fespace = mesh->GetNodes()->FESpace();
   }
   else
   {
      fec = new H1_FECollection(order, dim);
      fespace = new FiniteElementSpace(mesh, fec, dim);
   }


   if (!strcmp(geometry, "plate"))
   {
      fec_ho = new H1_FECollection(higher_order, dim);
      fespace_ho = new FiniteElementSpace(mesh, fec_ho, dim);
   }



   cout << "Number of finite element unknowns: " << fespace->GetTrueVSize()
        << endl << "Assembling: " << flush;

   // 6. Determine the list of true (i.e. conforming) essential boundary dofs.
   //    In this example, the boundary conditions are defined by marking only
   //    boundary attribute 1 from the mesh as essential and converting it to a
   //    list of true dofs.
   Array<int> ess_tdof_plist, ess_tdof_list, ess_bdr(mesh->bdr_attributes.Max());
   if (!strcmp(geometry, "beam"))
   {
      ess_bdr = 0;
      ess_bdr[20] = 1;
      fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }
   else if (!strcmp(geometry, "plate"))
   {
      ess_bdr = 0;
      ess_bdr[8] = 1;
      fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_plist, 0);
      ess_tdof_list.Append(ess_tdof_plist);
      ess_bdr = 0;
      ess_bdr[5] = 1;
      fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_plist, 1);
      ess_tdof_list.Append(ess_tdof_plist);
   }
   else if (!strcmp(geometry, "seat"))
   {
      ess_bdr = 0;
      ess_bdr[230] = 1;
      fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }

   // 7. Set up the linear form b(.) which corresponds to the right-hand side of
   //    the FEM linear system. In this case, b_i equals the boundary integral
   //    of f*phi_i where f represents a "pull down" force on the Neumann part
   //    of the boundary and phi_i are the basis functions in the finite element
   //    fespace. The force is defined by the VectorArrayCoefficient object f,
   //    which is a vector of Coefficient objects. The fact that f is non-zero
   //    on boundary attribute 2 is indicated by the use of piece-wise constants
   //    coefficient for its last component.
   VectorArrayCoefficient f(dim);
   for (int i = 0; i < dim-1; i++)
   {
      f.Set(i, new ConstantCoefficient(0.0));
   }
   {
      Vector pull_force(mesh->bdr_attributes.Max());
      pull_force = 0.0;
      if (!strcmp(geometry, "beam"))
      {
         pull_force(21) = -.2e-2;
         f.Set(dim-1, new PWConstCoefficient(pull_force));
      }
      else if (!strcmp(geometry, "plate"))
      {
         pull_force(6) = -10;
         f.Set(dim-2, new PWConstCoefficient(pull_force));
      }
      else if (!strcmp(geometry, "seat"))
      {
         pull_force(231) = -1;
         f.Set(dim-3, new PWConstCoefficient(pull_force));
      }
   }

   LinearForm *b = new LinearForm(fespace);
   b->AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(f));
   cout << "r.h.s. ... " << flush;
   b->Assemble();

   // 8. Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   GridFunction x(fespace);
   x = 0.0;

   // 9. Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the linear elasticity integrator with piece-wise
   //    constants coefficient lambda and mu.
   Vector lambda(mesh->attributes.Max());
   Vector mu(mesh->attributes.Max());
   if (!strcmp(geometry, "beam"))
   {
      lambda = 1.0;
      mu = 1.0;
   }
   else if (!strcmp(geometry, "plate"))
   {
      double E = 2e11, v = 0.3;
      lambda = (E * v)/((1+v)*(1-2.0*v));
      mu = E/(2.0*(1.0+v));
   }
   else if (!strcmp(geometry, "seat"))
   {
      lambda = 1.0;
      mu = 1.0;
   }
   PWConstCoefficient lambda_func(lambda);
   PWConstCoefficient mu_func(mu);

   BilinearForm *a = new BilinearForm(fespace);
   a->AddDomainIntegrator(new ElasticityIntegrator(lambda_func,mu_func));

   // 10. Assemble the bilinear form and the corresponding linear system,
   //     applying any necessary transformations such as: eliminating boundary
   //     conditions, applying conforming constraints for non-conforming AMR,
   //     static condensation, etc.
   cout << "matrix ... " << flush;
   if (static_cond) { a->EnableStaticCondensation(); }
   a->Assemble();

   SparseMatrix A;
   Vector B, X;
   a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);
   cout << "done." << endl;

   cout << "Size of linear system: " << A.Height() << endl;

#ifndef MFEM_USE_SUITESPARSE
   // 11. Define a simple symmetric Gauss-Seidel preconditioner and use it to
   //     solve the system Ax=b with PCG.
   GSSmoother M(A);
   PCG(A, M, B, X, 1, 500, 1e-8, 0.0);
#else
   // 11. If MFEM was compiled with SuiteSparse, use UMFPACK to solve the system.
   UMFPackSolver umf_solver;
   umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
   umf_solver.SetOperator(A);
   umf_solver.Mult(B, X);
#endif

   // 12. Recover the solution as a finite element grid function.
   a->RecoverFEMSolution(X, *b, x);

   if (!strcmp(geometry, "beam"))
   {
      VectorFunctionCoefficient D(dim, D_exact);
      double error = x.ComputeL2Error(D);
   
      cout << "\n|| d_h - d ||_{L^2} = " << error << '\n' << endl;
   }
   else if (!strcmp(geometry, "plate"))
   {
      GridFunction prolonged_coarse_function(fespace_ho);
      Operator* referenceOperator = nullptr;
      referenceOperator = new PRefinementTransferOperator(*fespace, *fespace_ho);
      referenceOperator->Mult(x, prolonged_coarse_function);

      ifstream mat_stream(solution_file);
      GridFunction ho_gf(mesh, mat_stream);
      ho_gf *= -1;
      VectorGridFunctionCoefficient ho_gf_co(&ho_gf);
      double error = prolonged_coarse_function.ComputeL2Error(ho_gf_co);
      cout << "\n|| d_h - d ||_{L^2} = " << error << '\n' << endl;
   }
   else if (!strcmp(geometry, "seat"))
   {
   }

   // 13. For non-NURBS meshes, make the mesh curved based on the finite element
   //     space. This means that we define the mesh elements through a fespace
   //     based transformation of the reference element. This allows us to save
   //     the displaced mesh as a curved mesh when using high-order finite
   //     element displacement field. We assume that the initial mesh (read from
   //     the file) is not higher order curved mesh compared to the chosen FE
   //     space.
   if (!mesh->NURBSext)
   {
      mesh->SetNodalFESpace(fespace);
   }

   // 14. Save the displaced mesh and the inverted solution (which gives the
   //     backward displacements to the original grid). This output can be
   //     viewed later using GLVis: "glvis -m displaced.mesh -g sol.gf".
   {
      GridFunction *nodes = mesh->GetNodes();
      std::cout << "No. of nodes on the grid: " << mesh->GetNodes()->Size() << std::endl;
      *nodes += x;
      x *= -1;
    
      std::string s = file;
      string meshfile_name = s + "mesh.";
      string solfile_name = s + "sol.";
       
      ostringstream mesh_name, sol_name;
      mesh_name << meshfile_name << setfill('0') << setw(6);
      sol_name << solfile_name << setfill('0') << setw(6);

      ofstream mesh_ofs(mesh_name.str().c_str());
      mesh_ofs.precision(8);
      mesh->Print(mesh_ofs);

      ofstream sol_ofs(sol_name.str().c_str());
      sol_ofs.precision(32);
      x.Save(sol_ofs);
   }

   // 15. Send the above data by socket to a GLVis server. Use the "n" and "b"
   //     keys in GLVis to visualize the displacements.
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock.precision(8);
      sol_sock << "solution\n" << *mesh << x << flush;
   }

   // 16. Free the used memory.
   delete a;
   delete b;
   if (fec)
   {
      delete fespace;
      delete fec;
   }
   delete mesh;

   return 0;
}

void D_exact(const Vector &x, Vector &D)
{
  double force, I, Dmod, lamda, mu, L;
  force = 0.002;
  lamda = 1.0;
  mu = 1.0;
  I = 1.0/12; 
  L = 8.0;
  Dmod = (mu * (3.0 * lamda + 2.0 * mu))/(lamda + mu);
    
  D(0) = 0.0;
  D(1) = 0.0;
  D(2) = -(force * x(0) * x(0)/(6.0 * Dmod * I)) * (3.0 * L - x(0));
}
