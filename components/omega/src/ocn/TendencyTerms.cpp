//===-- ocn/TendencyTerms.cpp - Tendency Terms ------------------*- C++ -*-===//
//
// The tendency terms that update state variables are implemented as functors,
// i.e. as classes that act like functions. This source defines the class
// constructors for these functors, which initialize the functor objects using
// the Mesh objects and info from the Config. The function call operators () are
// defined in the corresponding header file.
//
//===----------------------------------------------------------------------===//

#include "TendencyTerms.h"
#include "AuxiliaryState.h"
#include "DataTypes.h"
#include "HorzMesh.h"
#include "OceanState.h"
#include "Tracers.h"
#include "VertCoord.h"
#include "HorzOperators.h"

namespace OMEGA {

ThicknessFluxDivOnCell::ThicknessFluxDivOnCell(const HorzMesh *Mesh)
    : NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      DvEdge(Mesh->DvEdge), AreaCell(Mesh->AreaCell),
      EdgeSignOnCell(Mesh->EdgeSignOnCell) {}

PotentialVortHAdvOnEdge::PotentialVortHAdvOnEdge(const HorzMesh *Mesh)
    : NEdgesOnEdge(Mesh->NEdgesOnEdge), EdgesOnEdge(Mesh->EdgesOnEdge),
      WeightsOnEdge(Mesh->WeightsOnEdge), EdgeMask(Mesh->EdgeMask) {}

KEGradOnEdge::KEGradOnEdge(const HorzMesh *Mesh)
    : CellsOnEdge(Mesh->CellsOnEdge), DcEdge(Mesh->DcEdge),
      EdgeMask(Mesh->EdgeMask) {}

SSHGradOnEdge::SSHGradOnEdge(const HorzMesh *Mesh)
    : CellsOnEdge(Mesh->CellsOnEdge), DcEdge(Mesh->DcEdge),
      EdgeMask(Mesh->EdgeMask) {}

VelocityDiffusionOnEdge::VelocityDiffusionOnEdge(const HorzMesh *Mesh)
    : CellsOnEdge(Mesh->CellsOnEdge), VerticesOnEdge(Mesh->VerticesOnEdge),
      DcEdge(Mesh->DcEdge), DvEdge(Mesh->DvEdge),
      MeshScalingDel2(Mesh->MeshScalingDel2), EdgeMask(Mesh->EdgeMask) {}

VelocityHyperDiffOnEdge::VelocityHyperDiffOnEdge(const HorzMesh *Mesh)
    : CellsOnEdge(Mesh->CellsOnEdge), VerticesOnEdge(Mesh->VerticesOnEdge),
      DcEdge(Mesh->DcEdge), DvEdge(Mesh->DvEdge),
      MeshScalingDel4(Mesh->MeshScalingDel4), EdgeMask(Mesh->EdgeMask) {}

WindForcingOnEdge::WindForcingOnEdge(const HorzMesh *Mesh)
    : Enabled(false), EdgeMask(Mesh->EdgeMask) {}

BottomDragOnEdge::BottomDragOnEdge(const HorzMesh *Mesh,
                                   const VertCoord *VCoord)
    : Enabled(false), Coeff(0), CellsOnEdge(Mesh->CellsOnEdge),
      NVertLayers(VCoord->NVertLayers), EdgeMask(Mesh->EdgeMask) {}

TracerHorzAdvOnCell::TracerHorzAdvOnCell(const HorzMesh *Mesh)
    : NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      CellsOnEdge(Mesh->CellsOnEdge), EdgeSignOnCell(Mesh->EdgeSignOnCell),
      DvEdge(Mesh->DvEdge), AreaCell(Mesh->AreaCell), EdgeMask(Mesh->EdgeMask) {
}

TracerHighOrderHorzAdvOnCell::TracerHighOrderHorzAdvOnCell(const HorzMesh *Mesh)
    : HorzontalMesh(Mesh),
      NAdvCellsForEdge("NumberOfCellsContribToAdvectionAtEdge",
                       Mesh->NEdgesOwned),
      AdvCellsForEdge("IndexOfCellsContributingToAdvection", 
		      Mesh->NEdgesOwned, Mesh->MaxEdges2 + 2),
      AdvMaskHighOrder("MaskForHighOrderAdvectionTerms", Mesh->NEdgesAll),
      AdvCoefs("CommonAdvectionCoefficients", Mesh->MaxEdges2+2, Mesh->NEdgesAll),
      AdvCoefs3rd("CommonAdvectionCoeffsForHighOrder", Mesh->MaxEdges2+2,
                  Mesh->NEdgesAll),
      HighOrderFlxHorz("HigherOrderHorizontalFlux",Tracers::getNumTracers(), Mesh->NEdgesAll,Mesh->NVertLevels/VecLength),
      NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      CellsOnEdge(Mesh->CellsOnEdge), EdgeSignOnCell(Mesh->EdgeSignOnCell),
      DvEdge(Mesh->DvEdge), AreaCell(Mesh->AreaCell) {
   deepCopy(HighOrderFlxHorz, 0);
std::cout<<__FILE__<<":"<<__LINE__<<" "<<Tracers::getNumTracers()<<" "<< Mesh->NEdgesAll<<" "<<Mesh->NVertLevels/VecLength<<std::endl;
}

TracerDiffOnCell::TracerDiffOnCell(const HorzMesh *Mesh)
    : NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      CellsOnEdge(Mesh->CellsOnEdge), EdgeSignOnCell(Mesh->EdgeSignOnCell),
      DvEdge(Mesh->DvEdge), DcEdge(Mesh->DcEdge), AreaCell(Mesh->AreaCell),
      MeshScalingDel2(Mesh->MeshScalingDel2), EdgeMask(Mesh->EdgeMask) {}

TracerHyperDiffOnCell::TracerHyperDiffOnCell(const HorzMesh *Mesh)
    : NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      CellsOnEdge(Mesh->CellsOnEdge), EdgeSignOnCell(Mesh->EdgeSignOnCell),
      DvEdge(Mesh->DvEdge), DcEdge(Mesh->DcEdge), AreaCell(Mesh->AreaCell),
      MeshScalingDel4(Mesh->MeshScalingDel4), EdgeMask(Mesh->EdgeMask) {}

void TracerHighOrderHorzAdvOnCell::init() {
   const HorzMesh *Mesh   = this->HorzontalMesh;
   const auto MaxEdges2   = Mesh->MaxEdges2;
   const auto NEdgesAll   = Mesh->NEdgesAll;
   const auto NCellsOwned = Mesh->NCellsOwned;
   const auto NEdgesOwned = Mesh->NEdgesOwned;
   // Allocate Kokkos arrays in member data

   SecondDerivativeOnCell secondDerivativeOnCell(Mesh);
   Array3DReal DerivTwo("DerivTwo", MaxEdges2+2, 2, NEdgesAll);
   parallelFor(
       {NCellsOwned},
       KOKKOS_LAMBDA(int ICell) { secondDerivativeOnCell(DerivTwo, ICell); });
   // Compute masks and coefficients
   Kokkos::fence();
   MasksAndCoefficients masksAndCoefficients(Mesh, DerivTwo, NAdvCellsForEdge,
                                             AdvCellsForEdge, AdvMaskHighOrder,
                                             AdvCoefs, AdvCoefs3rd);
   Kokkos::fence();
   parallelFor(
       {NEdgesOwned},
       KOKKOS_LAMBDA(int IEdge) { masksAndCoefficients(IEdge); });
   Kokkos::fence();
}
} // end namespace OMEGA

//===----------------------------------------------------------------------===//
