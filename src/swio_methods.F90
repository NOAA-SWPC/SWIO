module swio_methods

  use ESMF
  use NUOPC
  use COMIO

  implicit none

  real(ESMF_KIND_R8), parameter :: earthRadius = 6.3712e+03_ESMF_KIND_R8

  private

  public :: SWIO_ArrayWrite
  public :: SWIO_FieldGetTimeStamp
  public :: SWIO_GridCreateLatLon
  public :: SWIO_MeshWriteCoord


  interface SWIO_ArrayWrite
    module procedure ArrayWrite
  end interface

  interface SWIO_FieldGetTimeStamp
    module procedure FieldGetTimeStamp
  end interface

  interface SWIO_GridCreateLatLon
    module procedure GridCreateLatLon
  end interface

  interface SWIO_MeshWriteCoord
    module procedure MeshWriteCoord
  end interface

contains

  subroutine CoordinateSetLinear(numValues, valueRange, values, minIndex, maxIndex, rc)
    integer,            intent(in)  :: numValues
    real(ESMF_KIND_R8), intent(in)  :: valueRange(:)
    real(ESMF_KIND_R8), pointer     :: values(:)
    integer, optional,  intent(in)  :: minIndex
    integer, optional,  intent(in)  :: maxIndex
    integer, optional,  intent(out) :: rc

    ! -- local variables
    integer            :: item, numItems
    integer            :: localMinIndex, localMaxIndex
    real(ESMF_KIND_R8) :: valueIncrement

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    if (.not.associated(values)) then
      call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
        msg="values pointer must be associated",&
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return
    end if

    values = 0._ESMF_KIND_R8

    if (numValues < 1) return

    localMinIndex = 1
    localMaxIndex = numValues
    if (present(minIndex)) localMinIndex = minIndex
    if (present(maxIndex)) localMaxIndex = maxIndex

    valueIncrement = 0._ESMF_KIND_R8

    if (numValues > 1) then
      numItems = numValues - 1
      valueIncrement = (valueRange(2) - valueRange(1)) / numItems
      do item = localMinIndex, localMaxIndex
        values(item) = valueRange(1) + (item - 1) * valueIncrement
      end do
      ! overwrite first and last element for accuracy
      if (localMinIndex ==         1) values(        1) = valueRange(1)
      if (localMinIndex == numValues) values(numValues) = valueRange(2)
    end if
    
  end subroutine CoordinateSetLinear

  subroutine ArrayWrite(array, io, name, tile, rc)
    type(ESMF_Array)               :: array
    class(COMIO_T)                 :: io
    character(len=*),  intent(in)  :: name
    integer, optional, intent(in)  :: tile
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer                              :: localrc, stat
    integer                              :: de, deCount, localDe, localDeCount
    integer                              :: i, item, lsize
    integer                              :: dimCount, rank, tileCount, tileId
    integer, dimension(:),   allocatable :: undistLBound, undistUBound
    integer, dimension(:),   allocatable :: distgridToArrayMap
    integer, dimension(:),   allocatable :: localDeToDeMap
    integer, dimension(:),   allocatable :: deToTileMap
    integer, dimension(:),   allocatable :: globalElemCount
    integer, dimension(:),   allocatable :: localElemCount, localElemStart
    integer, dimension(:,:), allocatable :: minIndexPTile, maxIndexPTile
    integer, dimension(:,:), allocatable :: minIndexPDe, maxIndexPDe
    real(ESMF_KIND_R8),          pointer :: fptr1d(:), fptr2d(:,:), fptr3d(:,:,:)
    type(ESMF_DistGrid)                  :: distgrid

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    tileId = 1
    if (present(tile)) tileId = tile

    call ESMF_ArrayGet(array, rank=rank, distgrid=distgrid, &
      dimCount=dimCount, tileCount=tileCount, deCount=deCount, &
      localDeCount=localDeCount, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    if (rank < 1 .or. rank > 3) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_OBJ_BAD, &
        msg="Only Arrays with rank=1,2,3 are supported", &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)
      return  ! bail out
    end if

    if (tileId < 1 .or. tileId > tileCount) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_OUTOFRANGE, &
        msg="tile number not found", &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)
      return  ! bail out
    end if

    allocate(minIndexPDe(dimCount,deCount), maxIndexPDe(dimCount,deCount),    &
      minIndexPTile(dimCount, tileCount), maxIndexPTile(dimCount, tileCount), &
      distgridToArrayMap(dimCount), localDeToDeMap(localDeCount), &
      deToTileMap(deCount), globalElemCount(rank), localElemCount(rank), &
      localElemStart(rank), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, msg="Unable to allocate memory", &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    ! retrieve minimum and maximum indices per DE
    call ESMF_DistgridGet(distgrid, minIndexPDe=minIndexPDe, &
      maxIndexPDe=maxIndexPDe, minIndexPTile=minIndexPTile, &
      maxIndexPTile=maxIndexPTile, localDeToDeMap=localDeToDeMap, &
      deToTileMap=deToTileMap, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    ! fill undistributed dimensions first
    if (rank > dimCount) then
      lsize = rank - dimCount
      allocate(undistLBound(lsize), undistUBound(lsize), stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, msg="Unable to allocate memory", &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out
      call ESMF_ArrayGet(array, undistLBound=undistLBound, undistUBound=undistUBound, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out

      globalElemCount     (dimCount+1:rank) = undistUBound - undistLBound + 1
      localElemStart  (dimCount+1:rank) = undistLBound
      localElemCount(dimCount+1:rank) = globalElemCount(dimCount+1:rank)
   
      deallocate(undistLBound, undistUBound, stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, msg="Unable to free memory", &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out
    end if
 
    ! now work on distributed dimensions
    call ESMF_ArrayGet(array, distgridToArrayMap=distgridToArrayMap, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    do localDe = 0, localDeCount - 1

      de = localDeToDeMap(localDe + 1) + 1

      if (tileId == deToTileMap(de)) then

        globalElemCount      = 0
        localElemCount = 0
        localElemStart   = 0
        do item = 1, dimCount
          i = distgridToArrayMap(item)
          if (i /= 0) then
            globalElemCount     (i) = maxIndexPTile(item,tileId) - minIndexPTile(item,tileId) + 1
            localElemStart  (i) = minIndexPDe(item,de)
            localElemCount(i) = maxIndexPDe(item,de) - minIndexPDe(item,de) + 1
          end if
        end do
   
        ! set dataset global and local domain
        call io % domain(globalElemCount, localElemStart, localElemCount)
        if (io % err % check(msg="Failure setting dataset domain", &
          line=__LINE__,  &
          file=__FILE__)) then
          call ESMF_LogSetError(ESMF_RC_FILE_UNEXPECTED, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
          return  ! bail out
        end if
    
        ! write dataset
        select case (rank)
          case (1)
            nullify(fptr1d)
            call ESMF_ArrayGet(array, localDe=localDe, farrayPtr=fptr1d, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__,  &
              file=__FILE__,  &
              rcToReturn=rc)) return  ! bail out
            call io % write(trim(name), fptr1d)
            if (io % err % check(msg="Failure writing dataset: "//trim(name), &
              line=__LINE__,  &
              file=__FILE__)) then
              call ESMF_LogSetError(ESMF_RC_FILE_WRITE, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, &
                file=__FILE__, &
                rcToReturn=rc)
              return  ! bail out
            end if
          case (2)
            nullify(fptr2d)
            call ESMF_ArrayGet(array, localDe=localDe, farrayPtr=fptr2d, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__,  &
              file=__FILE__,  &
              rcToReturn=rc)) return  ! bail out
            call io % write(trim(name), fptr2d)
            if (io % err % check(msg="Failure writing dataset: "//trim(name), &
              line=__LINE__,  &
              file=__FILE__)) then
              call ESMF_LogSetError(ESMF_RC_FILE_WRITE, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, &
                file=__FILE__, &
                rcToReturn=rc)
              return  ! bail out
            end if
          case (3)
            nullify(fptr3d)
            call ESMF_ArrayGet(array, localDe=localDe, farrayPtr=fptr3d, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__,  &
              file=__FILE__,  &
              rcToReturn=rc)) return  ! bail out
            call io % write(trim(name), fptr3d)
            if (io % err % check(msg="Failure writing dataset: "//trim(name), &
              line=__LINE__,  &
              file=__FILE__)) then
              call ESMF_LogSetError(ESMF_RC_FILE_WRITE, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, &
                file=__FILE__, &
                rcToReturn=rc)
              return  ! bail out
            end if
        end select
      end if

    end do

    deallocate(minIndexPDe, maxIndexPDe, minIndexPTile, maxIndexPTile, &
      distgridToArrayMap, localDeToDeMap, deToTileMap, globalElemCount,   &
      localElemCount, localElemStart, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, msg="Unable to free memory", &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

  end subroutine ArrayWrite

  subroutine MeshWriteCoord(mesh, io, meshloc, rc)
    type(ESMF_Mesh)                           :: mesh
    class(COMIO_T)                            :: io
    type(ESMF_MeshLoc), optional, intent(in)  :: meshloc
    integer,            optional, intent(out) :: rc

    ! local variables
    logical                                         :: isPresent
    integer                                         :: localrc, stat
    integer                                         :: item
    integer                                         :: de, deCount, dimCount, tileCount
    integer                                         :: numOwnedElements, numOwnedNodes
    integer                                         :: spatialDim
    integer,            dimension(1)                :: globalElemCount
    integer,            dimension(1)                :: localElemCount, localElemStart
    integer,            dimension(:,:), allocatable :: minIndexPDe,   maxIndexPDe
    integer,            dimension(:,:), allocatable :: minIndexPTile, maxIndexPTile
    real(ESMF_KIND_R8), dimension(:),   allocatable :: ownedCoords
    type(ESMF_DistGrid)                             :: distgrid
    type(ESMF_MeshLoc)                              :: lmeshloc
    
    integer :: localDeCount
    integer, dimension(:), allocatable :: localDeToDeMap

    character(len=*), parameter :: coordLabels(3) = &
      (/ "mesh_x", "mesh_y", "mesh_z" /)

    ! begin
    if (present(rc)) rc = ESMF_SUCCESS

    lmeshloc = ESMF_MESHLOC_NODE
    if (present(meshloc)) lmeshloc = meshloc

    if      (lmeshloc == ESMF_MESHLOC_NODE) then
      call ESMF_MeshGet(mesh, nodalDistgridIsPresent=isPresent, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out
      if (isPresent) then
        call ESMF_MeshGet(mesh, nodalDistgrid=distgrid, &
          numOwnedNodes=numOwnedNodes, spatialDim=spatialDim, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) return  ! bail out
        allocate(ownedCoords(spatialDim*numOwnedNodes), stat=stat)
        if (ESMF_LogFoundAllocError(statusToCheck=stat, &
          msg="Unable to allocate memory", &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) return  ! bail out
        call ESMF_MeshGet(mesh, ownedNodeCoords=ownedCoords, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) return  ! bail out
      else
        call ESMF_LogSetError(ESMF_RC_ARG_INCOMP, msg="nodal DistGrid not present", &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)
        return
      end if
    else if (lmeshloc == ESMF_MESHLOC_ELEMENT) then
      call ESMF_MeshGet(mesh, elementDistgridIsPresent=isPresent, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out
      if (isPresent) then
        call ESMF_MeshGet(mesh, elementDistgrid=distgrid, &
          numOwnedElements=numOwnedElements, spatialDim=spatialDim, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) return  ! bail out
        allocate(ownedCoords(spatialDim*numOwnedElements), stat=stat)
        if (ESMF_LogFoundAllocError(statusToCheck=stat, &
          msg="Unable to allocate memory", &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) return  ! bail out
        call ESMF_MeshGet(mesh, ownedElemCoords=ownedCoords, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) return  ! bail out
      else
        call ESMF_LogSetError(ESMF_RC_ARG_INCOMP, msg="element DistGrid not present", &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)
        return
      end if
    end if

    call ESMF_DistgridGet(distgrid, dimCount=dimCount, tileCount=tileCount, &
      deCount=deCount, localDeCount=localDeCount, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    if (localDeCount /= 1) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_OBJ_BAD, &
        msg="Expected localDeCount = 1 for Mesh objects", &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)
      return  ! bail out
    end if

    if (tileCount /= 1) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_OBJ_BAD, &
        msg="Expected tileCount = 1 for Mesh objects", &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)
      return  ! bail out
    end if

    allocate(minIndexPDe(dimCount, deCount), maxIndexPDe(dimCount, deCount),  &
      minIndexPTile(dimCount, tileCount), maxIndexPTile(dimCount, tileCount), &
      localDeToDeMap(localDeCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, msg="Unable to allocate memory", &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    ! retrieve minimum and maximum indices per DE
    call ESMF_DistgridGet(distgrid, minIndexPDe=minIndexPDe, &
      maxIndexPDe=maxIndexPDe, minIndexPTile=minIndexPTile, &
      maxIndexPTile=maxIndexPTile, localDeToDeMap=localDeToDeMap, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    ! now work on distributed dimensions
    de = localDeToDeMap(1) + 1
    globalElemCount(1) = maxIndexPTile(1,1) - minIndexPTile(1,1) + 1
    localElemStart (1) = minIndexPDe(1,de)
    localElemCount (1) = maxIndexPDe(1,de) - minIndexPDe(1,de) + 1
    
    ! set dataset global and local domain
    call io % domain(globalElemCount, localElemStart, localElemCount)
    if (io % err % check(msg="Failure setting dataset domain", &
      line=__LINE__,  &
      file=__FILE__)) then
      call ESMF_LogSetError(ESMF_RC_FILE_UNEXPECTED, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)
      return  ! bail out
    end if
    
    ! write dataset
    do item = 1, spatialDim
      call io % write(coordLabels(item), ownedCoords(item::spatialDim))
      if (io % err % check(msg="Failure writing dataset: "//coordLabels(item), &
        line=__LINE__,  &
        file=__FILE__)) then
        call ESMF_LogSetError(ESMF_RC_FILE_WRITE, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)
        return  ! bail out
      end if
    end do

    deallocate(minIndexPDe, maxIndexPDe, minIndexPTile, maxIndexPTile, &
      localDeToDeMap, ownedCoords, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, msg="Unable to free memory", &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

  end subroutine MeshWriteCoord

  subroutine FieldGetTimeStamp(field, timeStamp, isTimeValid, rc)
    type(ESMF_Field)               :: field
    character(len=15), intent(out) :: timeStamp
    logical,           intent(out) :: isTimeValid
    integer, optional, intent(out) :: rc

    ! local variables
    integer               :: localrc
    integer               :: mm, dd
    integer(ESMF_KIND_I4) :: yy, h, m, s
    type(ESMF_Time)       :: time

    ! begin
    if (present(rc)) rc = ESMF_SUCCESS

    timeStamp = ""

    ! get field time
    call NUOPC_GetTimestamp(field, isValid=isTimeValid, time=time, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    if (isTimeValid) then
      call ESMF_TimeGet(time, yy=yy, mm=mm, dd=dd, h=h, m=m, s=s, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      write(timeStamp, '(i4.4,2i2.2,"_",3i2.2)') yy, mm, dd, h, m, s
    end if

  end subroutine FieldGetTimeStamp

  function GridCreateLatLon(gcomp, rc) result (grid)
    type(ESMF_GridComp)            :: gcomp
    integer, optional, intent(out) :: rc
    type(ESMF_Grid)                :: grid

    ! local variables
    logical                         :: isLevelRelative
    integer                         :: localrc, stat
    integer                         :: dimCount, columnCount
    integer                         :: item
    integer                         :: verbosity
    integer                         :: lbnd(1), ubnd(1)
    integer, allocatable            :: dimLengths(:)
    character(len=ESMF_MAXSTR)      :: name
    character(len=ESMF_MAXSTR)      :: msgString
    real(ESMF_KIND_R8)              :: levelRange(2)
    real(ESMF_KIND_R8), pointer     :: coordPtr(:)
    real(ESMF_KIND_R8), allocatable :: levelValues(:)
    type(ESMF_Config)               :: config
    type(ESMF_StaggerLoc)           :: staggerloc

    ! local parameters
    character(len=*),   parameter :: rName = "GridCreateLatLon"
    real(ESMF_KIND_R8), parameter :: coordRange(2,2) = &
      reshape( (/   0._ESMF_KIND_R8, 360._ESMF_KIND_R8, &
                  -90._ESMF_KIND_R8,  90._ESMF_KIND_R8 /), &
               (/ 2, 2 /) )

    ! begin
    if (present(rc)) rc = ESMF_SUCCESS

    ! get component's info
    call NUOPC_CompGet(gcomp, name=name, verbosity=verbosity, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    ! get component's configuration
    call ESMF_GridCompGet(gcomp, config=config, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    ! retrieve grid size
    dimCount = ESMF_ConfigGetLen(config, label="output_grid_size:", rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return
    if (dimCount < 1 .or. dimCount > 3) then
      call ESMF_LogSetError(ESMF_RC_VAL_OUTOFRANGE, &
        msg="Invalid grid size: latlon supports 2 or 3 values", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return
    end if

    allocate(dimLengths(dimCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return

    ! get dimension lengths
    dimLengths = 0
    call ESMF_ConfigGetAttribute(config, dimLengths, &
      label="output_grid_size:", rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    if (dimCount == 3) then
      if (dimLengths(3) > 0) then
        call ESMF_ConfigGetAttribute(config, levelRange, &
          label="output_grid_level_range:", rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
      else
        ! get number of vertical levels provided in configuration file
        call ESMF_ConfigGetDim(config, dimLengths(3), columnCount, &
          label="output_grid_level_values::", rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
        allocate(levelValues(dimLengths(3)), stat=stat)
        if (ESMF_LogFoundAllocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return
        levelValues = 0._ESMF_KIND_R8
        ! read levels
        call ESMF_ConfigFindLabel(config, "output_grid_level_values::", rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
        do item = 1, dimLengths(3)
          call ESMF_ConfigNextLine(config, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) &
            return  ! bail out
          call ESMF_ConfigGetAttribute(config, levelValues(item), rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) &
            return  ! bail out
        end do
        ! verify that vertical coordinate is monotonically increasing
        if (any(levelValues(2:)-levelValues(1:dimLengths(3)-1) <= 0._ESMF_KIND_R8)) then
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
            msg="level values must increase strictly monotonically", &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)
          deallocate(levelValues, dimLengths, stat=stat)
          return
        end if
      end if
      ! check if vertical levels should be normalized
      call ESMF_ConfigGetAttribute(config, isLevelRelative, &
        label="output_grid_level_relative:", default=.false., rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
    end if

    ! create grid and add coordinates
    select case (dimCount)
      case (2)
        grid = ESMF_GridCreate1PeriDim( &
                 maxIndex  = dimLengths,&
                 coordSys  = ESMF_COORDSYS_SPH_DEG, &
                 coordDep1 = (/ 1 /), &
                 coordDep2 = (/ 2 /), &
                 indexFlag = ESMF_INDEX_GLOBAL, &
                 rc        = localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
        staggerloc = ESMF_STAGGERLOC_CENTER
      case (3)
        grid = ESMF_GridCreate1PeriDim( &
                 maxIndex  = dimLengths,&
                 coordSys  = ESMF_COORDSYS_SPH_DEG, &
                 coordDep1 = (/ 1 /), &
                 coordDep2 = (/ 2 /), &
                 coordDep3 = (/ 3 /), &
                 indexFlag = ESMF_INDEX_GLOBAL, &
                 rc        = localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
        staggerloc = ESMF_STAGGERLOC_CENTER_VCENTER
    end select

    ! add coordinates
    call ESMF_GridAddCoord(grid, staggerloc=staggerloc, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    ! set coordinate values
    do item = 1, dimCount
      nullify(coordPtr)
      call ESMF_GridGetCoord(grid, coordDim=item, staggerloc=staggerloc, &
        localDE=0, computationalLBound=lbnd, computationalUBound=ubnd, &
        farrayPtr=coordPtr, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      if (item == 3) then
        if (allocated(levelValues)) then
          coordPtr(lbnd(1):ubnd(1)) = levelValues(lbnd(1):ubnd(1))
        else
          call CoordinateSetLinear(dimLengths(item), levelRange, &
            coordPtr, minIndex=lbnd(1), maxIndex=ubnd(1), rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) &
            return  ! bail out
        end if
        if (isLevelRelative) then
          coordPtr(lbnd(1):ubnd(1)) = 1._ESMF_KIND_R8 + coordPtr(lbnd(1):ubnd(1)) / earthRadius
        end if
      else
        call CoordinateSetLinear(dimLengths(item), coordRange(:,item), &
          coordPtr, minIndex=lbnd(1), maxIndex=ubnd(1), rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
      end if
    end do

    if (btest(verbosity,8)) then
      write(msgString,'(a,": ",a,": lat/lon ",i0,"D grid created - resolution:",3(1x,i0))') &
        trim(name), rName, dimCount, dimLengths
      call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
    end if

    deallocate(dimLengths, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return

    if (allocated(levelValues)) then
      deallocate(levelValues, stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return
    end if

    ! store grid into component
    call ESMF_GridCompSet(gcomp, grid=grid, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

  end function GridCreateLatLon

end module swio_methods
