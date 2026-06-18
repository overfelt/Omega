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
#include "DataTypes.h"
#include "Error.h"
#include "HorzMesh.h"
#include "HorzOperators.h"
#include "Tracers.h"

namespace OMEGA {

PseudoThicknessFluxDivOnCell::PseudoThicknessFluxDivOnCell(
    const HorzMesh *Mesh, const VertCoord *VCoord)
    : NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      DvEdge(Mesh->DvEdge), AreaCell(Mesh->AreaCell),
      EdgeSignOnCell(Mesh->EdgeSignOnCell), MinLayerCell(VCoord->MinLayerCell),
      MaxLayerCell(VCoord->MaxLayerCell),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

PotentialVortHAdvOnEdge::PotentialVortHAdvOnEdge(const HorzMesh *Mesh,
                                                 const VertCoord *VCoord)
    : NEdgesOnEdge(Mesh->NEdgesOnEdge), EdgesOnEdge(Mesh->EdgesOnEdge),
      WeightsOnEdge(Mesh->WeightsOnEdge), EdgeMask(VCoord->EdgeMask),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

KEGradOnEdge::KEGradOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord)
    : CellsOnEdge(Mesh->CellsOnEdge), DcEdge(Mesh->DcEdge),
      EdgeMask(VCoord->EdgeMask), MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

SSHGradOnEdge::SSHGradOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord)
    : CellsOnEdge(Mesh->CellsOnEdge), DcEdge(Mesh->DcEdge),
      EdgeMask(VCoord->EdgeMask), MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

VelocityDiffusionOnEdge::VelocityDiffusionOnEdge(const HorzMesh *Mesh,
                                                 const VertCoord *VCoord)
    : CellsOnEdge(Mesh->CellsOnEdge), VerticesOnEdge(Mesh->VerticesOnEdge),
      DcEdge(Mesh->DcEdge), DvEdge(Mesh->DvEdge),
      MeshScalingDel2(Mesh->MeshScalingDel2), EdgeMask(VCoord->EdgeMask),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

VelocityHyperDiffOnEdge::VelocityHyperDiffOnEdge(const HorzMesh *Mesh,
                                                 const VertCoord *VCoord)
    : CellsOnEdge(Mesh->CellsOnEdge), VerticesOnEdge(Mesh->VerticesOnEdge),
      DcEdge(Mesh->DcEdge), DvEdge(Mesh->DvEdge),
      MeshScalingDel4(Mesh->MeshScalingDel4), EdgeMask(VCoord->EdgeMask),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

SfcStressForcingOnEdge::SfcStressForcingOnEdge(const HorzMesh *Mesh,
                                               const VertCoord *VCoord)
    : Enabled(false), EdgeMask(VCoord->EdgeMask),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot) {}

BottomDragOnEdge::BottomDragOnEdge(const HorzMesh *Mesh,
                                   const VertCoord *VCoord)
    : Enabled(false), Coeff(0), CellsOnEdge(Mesh->CellsOnEdge),
      NVertLayers(VCoord->NVertLayers), EdgeMask(VCoord->EdgeMask),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

TracerHorzAdvOnCell::TracerHorzAdvOnCell(const HorzMesh *Mesh,
                                         const VertCoord *VCoord,
					 const VertAdv *VAdv)
    : HorzontalMesh(Mesh), VerticalCoord(VCoord),
      NVertLayers(VCoord->NVertLayers),
      NAdvCellsForEdge("NumberOfCellsContribToAdvectionAtEdge",
                       Mesh->NEdgesAll),
      AdvCellsForEdge("IndexOfCellsContributingToAdvection", Mesh->NEdgesAll,
                      Mesh->MaxEdges2 + 2),
      AdvMaskHighOrder("MaskForHighOrderAdvectionTerms", Mesh->NEdgesAll,
                       VCoord->NVertLayers),
      CellsOnCell(Mesh->CellsOnCell),
      AdvCoefs("CommonAdvectionCoefficients", Mesh->MaxEdges2 + 2,
               Mesh->NEdgesAll),
      AdvCoefs3rd("CommonAdvectionCoeffsForHighOrder", Mesh->MaxEdges2 + 2,
                  Mesh->NEdgesAll),
      HighOrderFlxHorz("HigherOrderHorizontalFlux", Tracers::getNumTracers(),
                       Mesh->NEdgesAll, VCoord->NVertLayers),
      NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      CellsOnEdge(Mesh->CellsOnEdge), 
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MinLayerEdgeTop),
      EdgeSignOnCell(Mesh->EdgeSignOnCell),
      DvEdge(Mesh->DvEdge), AreaCell(Mesh->AreaCell),
      VerticalPseudoVelocity(VAdv->VerticalPseudoVelocity),
      HProvInv(),
      HNewInv(),
      HProv(),
      TracerMax(),
      TracerMin(),
      HighOrderFlx(),
      LowOrderFlx(),
      MinLayerCell(VCoord->MinLayerCell),
      MaxLayerCell(VCoord->MaxLayerCell),
      ActiveTracerHorizontalAdvectionEdgeFlux(),
      ActiveTracerHorizontalAdvectionTendency() {
   deepCopy(HighOrderFlxHorz, 0);
}

TracerDiffOnCell::TracerDiffOnCell(const HorzMesh *Mesh,
                                   const VertCoord *VCoord)
    : NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      CellsOnEdge(Mesh->CellsOnEdge), EdgeSignOnCell(Mesh->EdgeSignOnCell),
      DvEdge(Mesh->DvEdge), DcEdge(Mesh->DcEdge), AreaCell(Mesh->AreaCell),
      MeshScalingDel2(Mesh->MeshScalingDel2), EdgeMask(VCoord->EdgeMask),
      MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

TracerHyperDiffOnCell::TracerHyperDiffOnCell(const HorzMesh *Mesh,
                                             const VertCoord *VCoord)
    : NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      CellsOnEdge(Mesh->CellsOnEdge), EdgeSignOnCell(Mesh->EdgeSignOnCell),
      DvEdge(Mesh->DvEdge), DcEdge(Mesh->DcEdge), AreaCell(Mesh->AreaCell),
      MeshScalingDel4(Mesh->MeshScalingDel4), EdgeMask(VCoord->EdgeMask),
      MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

SurfaceTracerRestoringOnCell::SurfaceTracerRestoringOnCell(
    const HorzMesh *Mesh) {}

void TracerHorzAdvOnCell::init() {
   const HorzMesh *Mesh    = this->HorzontalMesh;
   const VertCoord *VCoord = this->VerticalCoord;
   const auto MaxEdges2    = Mesh->MaxEdges2;
   const auto NEdgesAll    = Mesh->NEdgesAll;
   const auto NCellsAll    = Mesh->NCellsAll;
   // Allocate Kokkos arrays in member data

   if (ForceLowOrder) {
      // Return when the 2nd-order tracer horz adv
      deepCopy(NAdvCellsForEdge, 0);
      deepCopy(AdvMaskHighOrder, 0);
      return;
   }

   SecondDerivativeOnCell secondDerivativeOnCell(Mesh);
   Array3DReal DerivTwo("DerivTwo", MaxEdges2 + 2, 2, NEdgesAll);
   parallelFor(
       {NCellsAll},
       KOKKOS_LAMBDA(int ICell) { secondDerivativeOnCell(DerivTwo, ICell); });
   // Compute masks and coefficients
   Kokkos::fence();
   MasksAndCoefficients masksAndCoefficients(
       Mesh, VCoord, DerivTwo, NAdvCellsForEdge, AdvCellsForEdge, 
       AdvMaskHighOrder, AdvCoefs, AdvCoefs3rd);
   Kokkos::fence();
   parallelFor(
       {NEdgesAll}, KOKKOS_LAMBDA(int IEdge) { masksAndCoefficients(IEdge); });
   Kokkos::fence();
   if (FCT) {
      HProvInv = Array2DReal("FCTProvesionalLayerThickness", Mesh->NEdgesAll,
                             NVertLayers);
      HNewInv =
          Array2DReal("FCTProvesionalNewInverse", Mesh->NEdgesAll, NVertLayers);
      HProv =
          Array2DReal("FCTProvesionalThickness", Mesh->NCellsAll, NVertLayers);
      TracerMax = Array2DReal("FCTTracerMax", Mesh->NCellsAll, NVertLayers);
      TracerMin = Array2DReal("FCTTracerMin", Mesh->NCellsAll, NVertLayers);
      HighOrderFlx = Array2DReal("FCTHighOrderFlx", Mesh->NEdgesAll, NVertLayers);
      LowOrderFlx = Array2DReal("FCTLowhOrderFlx", Mesh->NEdgesAll, NVertLayers);
      if (ComputeBudgets) {
	 const int NTracers          = Tracers::getNumTracers();
         ActiveTracerHorizontalAdvectionEdgeFlux = 
	   Array3DReal("FCTActiveTracerHorizontalAdvectionEdgeFlux", NTracers, Mesh->NCellsAll, NVertLayers);
         ActiveTracerHorizontalAdvectionTendency = 
	   Array3DReal("FCTActiveTracerHorizontalAdvectionTendency", NTracers, Mesh->NCellsAll, NVertLayers);
      }
   }
}
} // end namespace OMEGA

//===----------------------------------------------------------------------===//
