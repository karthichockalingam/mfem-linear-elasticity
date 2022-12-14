//                       MFEM Example 2 - Parallel Version
//
// Compile with: make ex2p
//
// Sample runs:  mpirun -np 4 ex2p -m ../data/beam-tri.mesh
//               mpirun -np 4 ex2p -m ../data/beam-quad.mesh
//               mpirun -np 4 ex2p -m ../data/beam-tet.mesh
//               mpirun -np 4 ex2p -m ../data/beam-hex.mesh
//               mpirun -np 4 ex2p -m ../data/beam-wedge.mesh
//               mpirun -np 4 ex2p -m ../data/beam-tri.mesh -o 2 -sys
//               mpirun -np 4 ex2p -m ../data/beam-quad.mesh -o 3 -elast
//               mpirun -np 4 ex2p -m ../data/beam-quad.mesh -o 3 -sc
//               mpirun -np 4 ex2p -m ../data/beam-quad-nurbs.mesh
//               mpirun -np 4 ex2p -m ../data/beam-hex-nurbs.mesh
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
   // 1. Initialize MPI and HYPRE.
   Mpi::Init(argc, argv);
   int num_procs = Mpi::WorldSize();
   int myid = Mpi::WorldRank();
   Hypre::Init();

   // 2. Parse command-line options.
   const char *mesh_file = "../data/beam-tri.mesh";
   const char *higher_order_mesh_file;
   int order = 1;
   int higher_order = 4;
   bool static_cond = false;
   bool visualization = 1;
   bool amg_elast = 0;
   bool reorder_space = false;
   const char *device_config = "cpu";
   int ref_levels = 0;
   const char *file = "xxx_";
   const char *solution_file;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&amg_elast, "-elast", "--amg-for-elasticity", "-sys",
                  "--amg-for-systems",
                  "Use the special AMG elasticity solver (GM/LN approaches), "
                  "or standard AMG for systems (unknown approach).");
   args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
                  "--no-static-condensation", "Enable static condensation.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&reorder_space, "-nodes", "--by-nodes", "-vdim", "--by-vdim",
                  "Use byNODES ordering of vector space instead of byVDIM");
   args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
   args.AddOption(&ref_levels, "-l", "--hlevel",
                  "Level of h refinement.");
   args.AddOption(&file, "-f", "--fileout", "Name of the output file.");
   args.AddOption(&solution_file, "-s2", "--solution2","Grid function for higher order solution.");
   args.AddOption(&higher_order_mesh_file, "-m2", "--mesh2",  "Mesh file for higher order solution.");
   args.AddOption(&higher_order, "-ho", "--higher_order","Finite element order (higher order polynomial degree).");
    
   args.Parse();
   if (!args.Good())
   {
      if (myid == 0)
      {
         args.PrintUsage(cout);
      }
      return 1;
   }
   if (myid == 0)
   {
      args.PrintOptions(cout);
   }

   // 3. Enable hardware devices such as GPUs, and programming models such as
   //    CUDA, OCCA, RAJA and OpenMP based on command line options.
   Device device(device_config);
   if (myid == 0) { device.Print(); }

   // 4. Read the (serial) mesh from the given mesh file on all processors.  We
   //    can handle triangular, quadrilateral, tetrahedral, hexahedral, surface
   //    and volume meshes with the same code.
   Mesh *mesh_higher_order = new Mesh(higher_order_mesh_file, 1, 1, false);
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();

   if (mesh->attributes.Max() < 2 || mesh->bdr_attributes.Max() < 2)
   {
      if (myid == 0)
         cerr << "\nInput mesh should have at least two materials and "
              << "two boundary attributes! (See schematic in ex2.cpp)\n"
              << endl;
      return 3;
   }

   // 5. Select the order of the finite element discretization space. For NURBS
   //    meshes, we increase the order by degree elevation.
   if (mesh->NURBSext)
   {
      mesh->DegreeElevate(order, order);
   }

   // 6. Refine the serial mesh on all processors to increase the resolution. In
   //    this example we do 'ref_levels' of uniform refinement. We choose
   //    'ref_levels' to be the largest number that gives a final mesh with no
   //    more than 1,000 elements.
   {
     cout << "Number of elements " << mesh->GetNE() << endl;
      // int ref_levels =
      //    (int)floor(log(1000./mesh->GetNE())/log(2.)/dim);
      cout << ">>> Refinement level: " << ref_levels << endl;
      for (int l = 0; l < ref_levels; l++)
      {
         mesh->UniformRefinement();
      }
   }

   // 7. Define a parallel mesh by a partitioning of the serial mesh. Refine
   //    this mesh further in parallel to increase the resolution. Once the
   //    parallel mesh is defined, the serial mesh can be deleted.
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   ParMesh *hmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;
   {
      int par_ref_levels = 0;
      for (int l = 0; l < par_ref_levels; l++)
      {
         pmesh->UniformRefinement();
      }
   }
   //ParMesh *hmesh = new ParMesh(*pmesh);
   // 8. Define a parallel finite element space on the parallel mesh. Here we
   //    use vector finite elements, i.e. dim copies of a scalar finite element
   //    space. We use the ordering by vector dimension (the last argument of
   //    the FiniteElementSpace constructor) which is expected in the systems
   //    version of BoomerAMG preconditioner. For NURBS meshes, we use the
   //    (degree elevated) NURBS space associated with the mesh nodes.
   FiniteElementCollection *fec, *higher_order_fec;
   ParFiniteElementSpace *fespace, *higher_order_fespace;

   const bool use_nodal_fespace = pmesh->NURBSext && !amg_elast;
   if (use_nodal_fespace)
   {
      fec = NULL;
      fespace = (ParFiniteElementSpace *)pmesh->GetNodes()->FESpace();
   }
   else
   {
      fec = new H1_FECollection(order, dim);
      higher_order_fec = new H1_FECollection(higher_order, dim);

      if (reorder_space)
      {
         std::cout << "Ordering byNodes " << std::endl;
         fespace = new ParFiniteElementSpace(pmesh, fec, dim, Ordering::byNODES);
         higher_order_fespace = new ParFiniteElementSpace(hmesh, higher_order_fec, dim, Ordering::byNODES);
      }
      else
      {
         std::cout << "Ordering byVDIM " << std::endl;
         fespace = new ParFiniteElementSpace(pmesh, fec, dim, Ordering::byVDIM);
         higher_order_fespace = new ParFiniteElementSpace(hmesh, higher_order_fec, dim, Ordering::byVDIM);
      }
   }
   HYPRE_BigInt size = fespace->GlobalTrueVSize();
   if (myid == 0)
   {
      cout << "Number of finite element unknowns: " << size << endl
           << "Assembling: " << flush;
   }

   // 9. Determine the list of true (i.e. parallel conforming) essential
   //    boundary dofs. In this example, the boundary conditions are defined by
   //    marking only boundary attribute 1 from the mesh as essential and
   //    converting it to a list of true dofs.
   cout << "# boundaries: " << pmesh->bdr_attributes.Max() << endl;
    Array<int> ess_tdof_plist, ess_tdof_list, ess_bdr(pmesh->bdr_attributes.Max());
    ess_bdr = 0;
    ess_bdr[8] = 1;
    fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_plist, 0);
    ess_tdof_list.Append(ess_tdof_plist);
    ess_bdr = 0;
    ess_bdr[5] = 1;
    fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_plist, 1);
    ess_tdof_list.Append(ess_tdof_plist);

   // 10. Set up the parallel linear form b(.) which corresponds to the
   //     right-hand side of the FEM linear system. In this case, b_i equals the
   //     boundary integral of f*phi_i where f represents a "pull down" force on
   //     the Neumann part of the boundary and phi_i are the basis functions in
   //     the finite element fespace. The force is defined by the object f, which
   //     is a vector of Coefficient objects. The fact that f is non-zero on
   //     boundary attribute 2 is indicated by the use of piece-wise constants
   //     coefficient for its last component.
   VectorArrayCoefficient f(dim);
   for (int i = 0; i < dim-1; i++)
   {
      f.Set(i, new ConstantCoefficient(0.0));
   }
   {
      Vector pull_force(pmesh->bdr_attributes.Max());
      pull_force = 0.0;
      pull_force(6) = -10;
      f.Set(0, new PWConstCoefficient(pull_force));
   }

   ParLinearForm *b = new ParLinearForm(fespace);
   b->AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(f));
   if (myid == 0)
   {
      cout << "r.h.s. ... " << flush;
   }
   b->Assemble();

   // 11. Define the solution vector x as a parallel finite element grid
   //     function corresponding to fespace. Initialize x with initial guess of
   //     zero, which satisfies the boundary conditions.
   ParGridFunction x(fespace);
   x = 0.0;

    // 12. Set up the bilinear form a(.,.) on the finite element space
    //    corresponding to the linear elasticity integrator with piece-wise
    //    constants coefficient lambda and mu.
    Vector lambda(pmesh->attributes.Max());
    double E = 2e11, v = 0.3;
    lambda = (E * v)/((1+v)*(1-2.0*v));
    //lambda(22) = lambda(23);
    PWConstCoefficient lambda_func(lambda);
    Vector mu(pmesh->attributes.Max());
    mu = E/(2.0*(1.0+v));
    //mu(22) = mu(23);
    PWConstCoefficient mu_func(mu);


   ParBilinearForm *a = new ParBilinearForm(fespace);
   a->AddDomainIntegrator(new ElasticityIntegrator(lambda_func, mu_func));

   // 13. Assemble the parallel bilinear form and the corresponding linear
   //     system, applying any necessary transformations such as: parallel
   //     assembly, eliminating boundary conditions, applying conforming
   //     constraints for non-conforming AMR, static condensation, etc.
   if (myid == 0) { cout << "matrix ... " << flush; }
   if (static_cond) { a->EnableStaticCondensation(); }
   a->Assemble();

   HypreParMatrix A;
   Vector B, X;
   a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);
   if (myid == 0)
   {
      cout << "done." << endl;
      cout << "Size of linear system: " << A.GetGlobalNumRows() << endl;
   }

   // 14. Define and apply a parallel PCG solver for A X = B with the BoomerAMG
   //     preconditioner from hypre.
   HypreBoomerAMG *amg = new HypreBoomerAMG(A);
   if (amg_elast && !a->StaticCondensationIsEnabled())
   {
      amg->SetElasticityOptions(fespace);
   }
   else
   {
      amg->SetSystemsOptions(dim, reorder_space);
   }
   HyprePCG *pcg = new HyprePCG(A);
   pcg->SetTol(1e-12);
   pcg->SetMaxIter(10000);
   pcg->SetPrintLevel(2);
   pcg->SetPreconditioner(*amg);
   pcg->Mult(B, X);

   // 15. Recover the parallel grid function corresponding to X. This is the
   //     local finite element solution on each processor.
   a->RecoverFEMSolution(X, *b, x);

   cout << "myid = " << myid << endl;

  
   PRefinementTransferOperator P(*fespace, *higher_order_fespace);
   GridFunction prolonged_coarse_function(higher_order_fespace);
   P.Mult(x, prolonged_coarse_function);
   ifstream mat_stream(solution_file);
   GridFunction solution_function(mesh_higher_order, mat_stream);
   GridFunctionCoefficient MeshSol(&solution_function);
   double error = prolonged_coarse_function.ComputeL2Error(MeshSol);


   if(myid == 0)
     cout << "\n|| d_h - d ||_{L^2} = " << error << '\n' << endl;
   
   
   // 16. For non-NURBS meshes, make the mesh curved based on the finite element
   //     space. This means that we define the mesh elements through a fespace
   //     based transformation of the reference element.  This allows us to save
   //     the displaced mesh as a curved mesh when using high-order finite
   //     element displacement field. We assume that the initial mesh (read from
   //     the file) is not higher order curved mesh compared to the chosen FE
   //     space.
   if (!use_nodal_fespace)
   {
      pmesh->SetNodalFESpace(fespace);
   }

   // 17. Save in parallel the displaced mesh and the inverted solution (which
   //     gives the backward displacements to the original grid). This output
   //     can be viewed later using GLVis: "glvis -np <np> -m mesh -g sol".
   {
      GridFunction *nodes = pmesh->GetNodes();
      *nodes += x;
      x *= -1;
    
      std::string s = file;
      string meshfile_name = s + "mesh.";
      string solfile_name = s + "sol.";
       
      ostringstream mesh_name, sol_name;
      mesh_name << meshfile_name << setfill('0') << setw(6) << myid;
      sol_name << solfile_name << setfill('0') << setw(6) << myid;

      ofstream mesh_ofs(mesh_name.str().c_str());
      mesh_ofs.precision(8);
      pmesh->Print(mesh_ofs);

      ofstream sol_ofs(sol_name.str().c_str());
      sol_ofs.precision(8);
      x.Save(sol_ofs);
   }

   // 18. Send the above data by socket to a GLVis server.  Use the "n" and "b"
   //     keys in GLVis to visualize the displacements.
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock << "parallel " << num_procs << " " << myid << "\n";
      sol_sock.precision(8);
      sol_sock << "solution\n" << *pmesh << x << flush;
   }

   // 19. Free the used memory.
   delete pcg;
   delete amg;
   delete a;
   delete b;
   if (fec)
   {
      delete fespace;
      delete fec;
   }
   delete pmesh;

   return 0;
}






void D_exact(const Vector &x, Vector &D)
{
    double T = 10.0;
    double a = 0.4;
    double E = 2e11;
    double v = 0.3;
    double mu = E/(2.0*(1.0+v));
    double r = std::sqrt(x(0)*x(0) + x(1)*x(1));
    double theta = std::atan2(x(1),x(0));
    double k = (3.0-v)/(1.0+v);

    double factor = T * a/(8.0 * mu);
    
    double term_one = (r/a)*(k+1.0)*std::cos(theta);
    double term_two = (2.0*a/r)*((k+1.0)*std::cos(theta)+std::cos(3.0*theta));
    double term_three = -(2.0*a*a*a/(r*r*r))*std::cos(3.0*theta);
    
    D(0) = factor * (term_one + term_two + term_three);
    
    term_one = (r/a)*(k-3.0)*std::sin(theta);
    term_two = (2.0*a/r)*((1.0-k)*std::sin(theta)+std::sin(3.0*theta));
    term_three = -(2.0*a*a*a/(r*r*r))*std::sin(3.0*theta);
    
    D(1) = factor * (term_one + term_two + term_three);
    D(2) = 0.0;
  
}

