module swio_methods

  use ESMF
  use NUOPC
  use COMIO

  use swio_data

  implicit none

  real(ESMF_KIND_R8), parameter :: earthRadius = 6.3712e+03_ESMF_KIND_R8

  private

  public :: SWIO_ArrayWrite
  public :: SWIO_FieldWrite
  public :: SWIO_FieldGetTimeStamp
  public :: SWIO_FieldWriteCoord
  public :: SWIO_GridCreateLatLon
  public :: SWIO_MeshWriteCoord
  public :: SWIO_MetadataParse
  public :: SWIO_OutputFieldsParse
  public :: SWIO_OutputMaskParse
  public :: SWIO_Output


  interface SWIO_ArrayWrite
    module procedure ArrayWrite
  end interface

  interface SWIO_FieldWrite
    module procedure FieldWrite
  end interface

  interface SWIO_FieldGetTimeStamp
    module procedure FieldGetTimeStamp
  end interface

  interface SWIO_FieldWriteCoord
    module procedure FieldWriteCoord
  end interface

  interface SWIO_GridCreateLatLon
    module procedure GridCreateLatLon
  end interface

  interface SWIO_MeshWriteCoord
    module procedure MeshWriteCoord
  end interface

  interface SWIO_MetadataParse
    module procedure MetadataParse
  end interface

  interface SWIO_OutputFieldsParse
    module procedure OutputFieldsParse
  end interface

  interface SWIO_OutputMaskParse
    module procedure OutputMaskParse
  end interface

  interface SWIO_Output
    module procedure FileWrite
  end interface

contains

  subroutine CoordinateSetLinear(numValues, valueRange, values, minIndex, maxIndex, wrap, rc)
    integer,            intent(in)  :: numValues
    real(ESMF_KIND_R8), intent(in)  :: valueRange(:)
    real(ESMF_KIND_R8), pointer     :: values(:)
    integer, optional,  intent(in)  :: minIndex
    integer, optional,  intent(in)  :: maxIndex
    logical, optional,  intent(in)  :: wrap
    integer, optional,  intent(out) :: rc

    ! -- local variables
    logical            :: addEndpoint
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

    addEndpoint = .true.
    if (present(wrap)) addEndpoint = .not.wrap

    values = 0._ESMF_KIND_R8

    if (numValues < 1) return

    localMinIndex = 1
    localMaxIndex = numValues
    if (present(minIndex)) localMinIndex = minIndex
    if (present(maxIndex)) localMaxIndex = maxIndex

    valueIncrement = 0._ESMF_KIND_R8

    if (numValues > 1) then
      numItems = numValues
      if (addEndpoint) numItems = numItems - 1
      valueIncrement = (valueRange(2) - valueRange(1)) / numItems
      do item = localMinIndex, localMaxIndex
        values(item) = valueRange(1) + (item - 1) * valueIncrement
      end do
      ! overwrite first and last element for accuracy
      if (localMinIndex == 1) values(1) = valueRange(1)
      if ((localMinIndex == numValues) .and. addEndpoint) &
        values(numValues) = valueRange(2)
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
    globalElemCount = 0
    localElemCount  = 0
    localElemStart  = 0
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

      globalElemCount(dimCount+1:rank) = undistUBound - undistLBound + 1
      localElemStart (dimCount+1:rank) = undistLBound
      localElemCount (dimCount+1:rank) = globalElemCount(dimCount+1:rank)
   
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

        do item = 1, dimCount
          i = distgridToArrayMap(item)
          if (i /= 0) then
            globalElemCount(i) = maxIndexPTile(item,tileId) - minIndexPTile(item,tileId) + 1
            localElemStart (i) = minIndexPDe(item,de)
            localElemCount (i) = maxIndexPDe(item,de) - minIndexPDe(item,de) + 1
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

  subroutine FieldWrite(field, io, tile, name, phaseName, rc)
    type(ESMF_Field)                        :: field
    class(COMIO_T)                          :: io
    integer,          optional, intent(in)  :: tile
    character(len=*), optional, intent(in)  :: name
    character(len=*), optional, intent(in)  :: phaseName
    integer,          optional, intent(out) :: rc

    ! local variables
    integer                    :: localrc
    integer                    :: item
    logical                    :: isSet
    character(len=ESMF_MAXSTR) :: pName
    character(len=ESMF_MAXSTR) :: fieldName
    character(len=ESMF_MAXSTR) :: svalue
    type(ESMF_Array)           :: array

    ! local parameters
    character(len=*), parameter :: attributeList(2,2) = &
      reshape( (/ &
        "LongName ", "long_name", &
        "Units    ", "units    "  &
      /), (/2,2/), order = (/2,1/) )

    ! begin
    if (present(rc)) rc = ESMF_SUCCESS

    pName = "FieldWrite"
    if (present(phaseName)) pName = phaseName

    ! get field info
     call ESMF_FieldGet(field, name=fieldName, array=array, rc=localrc)
     if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__,  &
       file=__FILE__,  &
       rcToReturn=rc)) &
       return  ! bail out

    ! replace output name if requested
    if (present(name)) fieldName = name

    call ArrayWrite(array, io, trim(fieldName), tile=tile, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    ! write attributes
    do item = 1, size(attributeList, dim=1)
      call NUOPC_GetAttribute(field, name=trim(attributeList(item,1)), &
        value=svalue, isSet=isSet, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      if (isSet) then
        call io % describe(trim(fieldName), trim(attributeList(item,2)), trim(svalue))
        if (io % err % check(msg="Failure writing attribute " &
          //trim(attributeList(item,2))//" for "//trim(fieldName), &
          line=__LINE__,  &
          file=__FILE__)) then
          call ESMF_LogSetError(ESMF_RC_FILE_WRITE, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
          return  ! bail out
        end if
      end if
    end do

  end subroutine FieldWrite

  subroutine MeshWriteCoord(mesh, io, meshloc, rc)
    type(ESMF_Mesh)                           :: mesh
    class(COMIO_T)                            :: io
    type(ESMF_MeshLoc), optional, intent(in)  :: meshloc
    integer,            optional, intent(out) :: rc

    ! local variables
    logical                                         :: isMemFreed, isPresent
    integer                                         :: localrc, stat
    integer                                         :: item, is, ie
    integer                                         :: de, deCount, localDe, localDeCount
    integer                                         :: dimCount, tileCount
    integer                                         :: numOwnedElements, numOwnedNodes
    integer                                         :: spatialDim
    integer,            dimension(1)                :: globalElemCount
    integer,            dimension(1)                :: localElemCount, localElemStart
    integer,            dimension(:),   allocatable :: localDeToDeMap
    integer,            dimension(:,:), allocatable :: minIndexPDe,   maxIndexPDe
    integer,            dimension(:,:), allocatable :: minIndexPTile, maxIndexPTile
    real(ESMF_KIND_R8), dimension(:),   allocatable :: ownedCoords
    type(ESMF_DistGrid)                             :: distgrid
    type(ESMF_MeshLoc)                              :: lmeshloc

    character(len=*), parameter :: coordLabels(3) = &
      (/ "mesh_x", "mesh_y", "mesh_z" /)

    ! begin
    if (present(rc)) rc = ESMF_SUCCESS

    ! check if coordinate and connection information is present
    call ESMF_MeshGet(mesh, isMemFreed=isMemFreed, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    ! return if not present
    if (isMemFreed) return

    ! set mesh location
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

    if (dimCount /= 1) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_OBJ_BAD, &
        msg="Expected dimCount = 1 for Mesh objects", &
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
    globalElemCount = 0
    localElemStart  = 0
    localElemCount  = 0

    is = 0
    do localDe = 0, localDeCount - 1

      de = localDeToDeMap(localDe + 1) + 1
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
      ie = is + spatialDim * localElemCount(1)
      do item = 1, spatialDim
        call io % write(coordLabels(item), ownedCoords(is+item:ie:spatialDim))
        if (io % err % check(msg="Failure writing dataset: "//coordLabels(item), &
          line=__LINE__,  &
          file=__FILE__)) then
          call ESMF_LogSetError(ESMF_RC_FILE_WRITE, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)
          return  ! bail out
        end if
      end do
      is = ie
    end do

    deallocate(minIndexPDe, maxIndexPDe, minIndexPTile, maxIndexPTile, &
      localDeToDeMap, ownedCoords, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, msg="Unable to free memory", &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

  end subroutine MeshWriteCoord

  subroutine FieldWriteCoord(field, io, georeference, logLabel, verbosity, rc)
    type(ESMF_Field)                        :: field
    class(COMIO_T)                          :: io
    logical,          optional, intent(in)  :: georeference
    character(len=*), optional, intent(in)  :: loglabel
    integer,          optional, intent(in)  :: verbosity
    integer,          optional, intent(out) :: rc

    ! local variables
    integer                    :: localrc, stat
    integer                    :: dimCount, rank
    integer                    :: iCoord, item, n
    integer                    :: lverbosity
    integer                    :: dimLength(1)
    integer, allocatable       :: ungriddedLBound(:), ungriddedUBound(:)
    logical                    :: lgeoreference
    character(len=ESMF_MAXSTR) :: llabel
    character(len=ESMF_MAXSTR) :: msgString
    type(ESMF_Array)           :: array
    type(ESMF_GeomType_Flag)   :: geomtype
    type(ESMF_Grid)            :: grid
    type(ESMF_Mesh)            :: mesh
    type(ESMF_MeshLoc)         :: meshloc

    ! local parameters
    character(len=*), parameter :: coordLabels(3,2) = &
      reshape( (/ &
        "grid_x", "grid_y", "grid_z", &
        "lon   ", "lat   ", "alt   "  &
      /), (/ 3, 2 /) )
    character(len=*), parameter :: geoAttributes(4,3) = &
      reshape( (/ &
        "long_name    ", "longitude    ", &
        "units        ", "degrees_east ", &
        "long_name    ", "latitude     ", &
        "units        ", "degrees_north", &
        "long_name    ", "altitude     ", &
        "units        ", "km           "  &
      /), (/ 4, 3 /) )

    ! begin
    if (present(rc)) rc = ESMF_SUCCESS

    ! set defaults
    llabel = "FieldWriteCoord"
    if (present(logLabel)) llabel = logLabel

    lverbosity = 0
    if (present(verbosity)) lverbosity = verbosity

    lgeoreference = .false.
    if (present(georeference)) lgeoreference = georeference

    ! get field GeomType
    call ESMF_FieldGet(field, geomtype=geomtype, rank=rank, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    if (geomtype == ESMF_GEOMTYPE_GRID) then
      ! write grid coordinates
      call ESMF_FieldGet(field, grid=grid, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      call ESMF_GridGet(grid, dimCount=dimCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out
      ! select set of coordinate labels
      iCoord = 1
      if (lgeoreference) iCoord = 2
      do item = 1, dimCount
        call ESMF_GridGetCoord(grid, item, array=array, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
        call SWIO_ArrayWrite(array, io, coordLabels(item,iCoord), rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
        if (lgeoreference) then
          do n = 1, size(geoAttributes, dim=1), 2
            call io % describe(trim(coordLabels(item,iCoord)), &
              trim(geoAttributes(n,item)), trim(geoAttributes(n + 1,item)))
            if (io % err % check(msg="Failure writing attribute " &
              //trim(geoAttributes(n,item))//" for "//trim(coordLabels(item,iCoord)), &
              line=__LINE__,  &
              file=__FILE__)) then
              call ESMF_LogSetError(ESMF_RC_FILE_WRITE, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__, &
                file=__FILE__, &
                rcToReturn=rc)
              return  ! bail out
            end if
          end do
        else
          call SWIO_ArrayWrite(array, io, coordLabels(item,iCoord), rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) &
            return  ! bail out
        end if
        if (btest(lverbosity,8)) then
          call ESMF_LogWrite(trim(llabel)//": Written coordinate "&
            //coordLabels(item,iCoord), ESMF_LOGMSG_INFO, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) &
            return  ! bail out
        end if
      end do
    else if (geomtype == ESMF_GEOMTYPE_MESH) then
      ! get mesh location the Field is built on
      call ESMF_FieldGet(field, mesh=mesh, meshloc=meshloc, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return
      ! set dimCount to 1 for Mesh objects
      dimCount = 1
      call SWIO_MeshWriteCoord(mesh, io, meshloc=meshloc, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      if (btest(lverbosity,8)) then
        call ESMF_LogWrite(trim(llabel)//": Written coordinates ",&
          ESMF_LOGMSG_INFO, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
      end if
    end if

    if (rank > dimCount) then
      n = rank - dimCount
      allocate(ungriddedLBound(n), ungriddedUBound(n), stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      call ESMF_FieldGet(field, ungriddedLBound=ungriddedLBound, &
        ungriddedUBound=ungriddedUBound, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      do item = 1, n
        dimLength(1) = ungriddedUBound(item) - ungriddedLBound(item) + 1
        call io % domain(dimLength, (/ 1 /), dimLength)
        if (io % err % check(msg="Failure setting dataset ungridded domain", &
          line=__LINE__,  &
          file=__FILE__)) then
          call ESMF_LogSetError(ESMF_RC_FILE_UNEXPECTED, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
          return  ! bail out
        end if
        if (btest(lverbosity,8)) then
          write(msgString,'(a,1x,i0)') trim(llabel)//": Written ungridded dimension", item
          call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) &
            return  ! bail out
        end if
      end do
      deallocate(ungriddedLBound, ungriddedUBound, stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
    end if

  end subroutine FieldWriteCoord

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

  subroutine FieldSetMask(field, maskField, maskValue, fillValue, rc)
    type(ESMF_Field)                :: field
    type(ESMF_Field),   intent(in)  :: maskField
    real(ESMF_KIND_R8), intent(in)  :: maskValue
    real(ESMF_KIND_R8), intent(in)  :: fillValue
    integer, optional,  intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: l, localDe, localDeCount, mrank, rank
    real(ESMF_KIND_R8), dimension(:,:),   pointer :: mptr2d
    real(ESMF_KIND_R8), dimension(:,:),   pointer :: fptr2d
    real(ESMF_KIND_R8), dimension(:,:,:), pointer :: fptr3d

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field, localDeCount=localDeCount, rank=rank, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    call ESMF_FieldGet(maskField, rank=mrank, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
       line=__LINE__,  &
       file=__FILE__,  &
       rcToReturn=rc)) &
       return  ! bail out

    do localDe = 0, localDeCount-1
      nullify(mptr2d)
      select case (mrank)
        case (2)
          call ESMF_FieldGet(maskField, localDe=localDe, farrayPtr=mptr2d, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) return  ! bail out
        case (3)
          nullify(fptr3d)
          call ESMF_FieldGet(maskField, localDe=localDe, farrayPtr=fptr3d, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) return  ! bail out
          mptr2d => fptr3d(:,:,lbound(fptr3d,dim=3))
      end select
      if (associated(mptr2d)) then
        select case (rank)
          case (2)
            nullify(fptr2d)
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=fptr2d, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__,  &
              file=__FILE__,  &
              rcToReturn=rc)) return  ! bail out
            where (mptr2d == maskValue) fptr2d = fillValue
          case (3)
            nullify(fptr3d)
            call ESMF_FieldGet(field, localDe=localDe, farrayPtr=fptr3d, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__,  &
              file=__FILE__,  &
              rcToReturn=rc)) return  ! bail out
            ! -- set mask from top to bottom to prevent issues when mask
            ! -- and masked fields are the same object
            do l = ubound(fptr3d,dim=3), lbound(fptr3d,dim=3), -1
              where (mptr2d == maskValue) fptr3d(:,:,l) = fillValue
            end do
        end select
      end if
    end do

  end subroutine FieldSetMask

  function GridCreateLatLon(gcomp, rc) result (grid)
    type(ESMF_GridComp)            :: gcomp
    integer, optional, intent(out) :: rc
    type(ESMF_Grid)                :: grid

    ! local variables
    logical                         :: isLevelRelative
    integer                         :: localrc, stat
    integer                         :: dimCount, columnCount
    integer                         :: ungriddedDimLength
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
    ungriddedDimLength = 0
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
      else if (dimLengths(3) == 0) then
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
      else
        ! 3rd dimension is ungridded
        dimCount = 2
        ungriddedDimLength = -dimLengths(3)
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
                 maxIndex  = dimLengths(1:dimCount),&
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
          coordPtr, minIndex=lbnd(1), maxIndex=ubnd(1), wrap=(item==1), rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
      end if
    end do

    ! add ungridded dimension if needed
    if (ungriddedDimLength > 0) then
      ! record ungridded dimension size as ESMF Attribute
      call ESMF_AttributeAdd(grid, convention="NUOPC", purpose="Instance", &
        attrList=(/ "UngriddedDimLength" /), rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

      call ESMF_AttributeSet(grid, name="UngriddedDimLength", &
        value=ungriddedDimLength, convention="NUOPC", purpose="Instance", rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
    end if

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

  subroutine FileWrite(gcomp, phaseName, rc)
    type(ESMF_GridComp)                     :: gcomp
    character(len=*), optional, intent(in)  :: phaseName
    integer,          optional, intent(out) :: rc

    ! local variables
    logical                             :: isTimeValid
    integer                             :: localrc
    integer                             :: i, item, stat
    integer                             :: dimCount
    integer                             :: fieldCount
    integer                             :: itemCount
    integer                             :: verbosity
    character(len=15)                   :: timeStamp
    character(len=ESMF_MAXSTR)          :: name
    character(len=ESMF_MAXSTR)          :: pName
    character(len=ESMF_MAXSTR)          :: fieldName
    character(len=ESMF_MAXSTR)          :: fileName
    character(len=ESMF_MAXSTR), pointer :: standardNameList(:)
    type(ESMF_Array)                    :: array
    type(ESMF_Field)                    :: field
    type(ESMF_Field)                    :: maskField
    type(ESMF_GeomType_Flag)            :: geomtype
    type(ESMF_Grid)                     :: grid
    type(ESMF_Mesh)                     :: mesh
    type(ESMF_MeshLoc)                  :: meshloc
    type(ESMF_State)                    :: importState
    type(SWIO_InternalState_T)          :: is
    type(SWIO_Data_T), pointer          :: this

    ! begin
    if (present(rc)) rc = ESMF_SUCCESS

    pName = "FileWrite"
    if (present(phaseName)) pName = phaseName

    ! get component's info
    call NUOPC_CompGet(gcomp, name=name, verbosity=verbosity, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    ! get component's internal state
    call ESMF_GridCompGetInternalState(gcomp, is, localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out
    this => is % wrap

    if (this % outputCount == 0) then
      if (btest(verbosity,8)) then
        call ESMF_LogWrite(trim(name)//": "//trim(pName)//": No fields to write", &
          ESMF_LOGMSG_INFO, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
      end if
      return
    end if

    ! query the component for its importState
    call ESMF_GridCompGet(gcomp, importState=importState, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    ! retrieve field list
    nullify(standardNameList)
    call NUOPC_GetStateMemberLists(importState, &
      StandardNameList=standardNameList, nestedFlag=.true., rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    itemCount = 0
    if (associated(standardNameList)) itemCount = size(standardNameList)

    if (itemCount == 0) then
      if (btest(verbosity,8)) then
        call ESMF_LogWrite(trim(name)//": "//trim(pName)//": No imported fields", &
          ESMF_LOGMSG_INFO, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
      end if
      return
    end if

    ! -- check if time stamp of imported fields is valid -- use field #1
    item = 1
    call ESMF_StateGet(importState, trim(standardNameList(item)), field, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    call SWIO_FieldGetTimeStamp(field, timeStamp, isTimeValid, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    if (isTimeValid) then
      ! open file
      fileName = trim(this % filePrefix) // "." // timeStamp // "." // trim(this % fileSuffix)
      call this % io % open(fileName, this % cmode)
      if (this % io % err % check(msg="Failure creating file "//trim(fileName), &
        line=__LINE__,  &
        file=__FILE__)) then
        call ESMF_LogSetError(ESMF_RC_FILE_CREATE, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__,  &
          rcToReturn=rc)
        return  ! bail out
      end if
      if (btest(verbosity,8)) then
        call ESMF_LogWrite(trim(name)//": "//trim(pName)//": Opened file "//trim(fileName), &
          ESMF_LOGMSG_INFO, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
      end if

      ! write Field coordinates and ungridded dimensions
      call FieldWriteCoord(field, this % io, georeference=this % geoReference, &
        logLabel=trim(name)//": "//trim(pName), verbosity=verbosity, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

      ! mask output fields using fill value
      if (associated(this % mask)) then
        itemCount = 0
        if (associated(this % output)) itemCount = size(this % output)
        do item = 1, itemCount
          call ESMF_StateGet(importState, trim(this % output(item) % key), &
            field, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) &
            return  ! bail out
          if (field /= this % mask % field) then
            call FieldSetMask(field, this % mask % field, this % mask % value, &
              this % mask % fill, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__,  &
              file=__FILE__,  &
              rcToReturn=rc)) &
              return  ! bail out
          end if
        end do
        itemCount = 0
        if (associated(this % task)) itemCount = size(this % task)
        do item = 1, itemCount
          fieldCount = 0
          if (associated(this % task(item) % fieldOut)) &
            fieldCount = size(this % task(item) % fieldOut)
          do i = 1, fieldCount
            call FieldSetMask(this % task(item) % fieldOut(i), this % mask % field, &
              this % mask % value, this % mask % fill, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__,  &
              file=__FILE__,  &
              rcToReturn=rc)) &
              return  ! bail out
          end do
        end do
        call FieldSetMask(this % mask % field, this % mask % field, &
          this % mask % value, this % mask % fill, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
        call this % io % fill(this % mask % fill)
        if (this % io % err % check(msg="Failure setting fill value in file " &
          //trim(fileName), &
          line=__LINE__,  &
          file=__FILE__)) then
          call ESMF_LogSetError(ESMF_RC_FILE_CREATE, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__,  &
            rcToReturn=rc)
           return  ! bail out
        end if
      end if

      ! write select output fields
      itemCount = 0
      if (associated(this % output)) itemCount = size(this % output)

      do item = 1, itemCount
        call ESMF_StateGet(importState, trim(this % output(item) % key), field, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
        call FieldWrite(field, this % io, phaseName=pName, &
          name=this % output(item) % value, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
        if (btest(verbosity,8)) then
          call ESMF_LogWrite(trim(name)//": "//trim(pName)//": Written "&
            //trim(standardNameList(item)), ESMF_LOGMSG_INFO, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) &
            return  ! bail out
        end if
      end do

      ! write computed fields if present
      itemCount = 0
      if (associated(this % task)) itemCount = size(this % task)
      do item = 1, itemCount
        fieldCount = 0
        if (associated(this % task(item) % fieldOut)) &
          fieldCount = size(this % task(item) % fieldOut)
        do i = 1, fieldCount
          call FieldWrite(this % task(item) % fieldOut(i), this % io, &
            phaseName=pName, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) &
            return  ! bail out
          if (btest(verbosity,8)) then
            call ESMF_FieldGet(this % task(item) % fieldOut(i), name=fieldName, &
              rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__,  &
              file=__FILE__,  &
              rcToReturn=rc)) &
              return  ! bail out
            call ESMF_LogWrite(trim(name)//": "//trim(pName)//": Written "&
              //trim(fieldName), ESMF_LOGMSG_INFO, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__,  &
              file=__FILE__,  &
              rcToReturn=rc)) &
              return  ! bail out
          end if
        end do
      end do

      ! write metadata as global attributes if provided
      itemCount = 0
      if (associated(this % meta)) itemCount = size(this % meta)
      do item = 1, itemCount
        if (trim(this % meta(item) % value) == "__swio_field_timestamp__") then
          call this % io % describe(trim(this % meta(item) % key), timeStamp)
        else
          call this % io % describe(trim(this % meta(item) % key), trim(this % meta(item) % value))
        end if
        if (this % io % err % check(msg="Failure writing global attribute " &
          //trim(this % meta(item) % key)//" to "//trim(fileName), &
          line=__LINE__,  &
          file=__FILE__)) then
          call ESMF_LogSetError(ESMF_RC_FILE_WRITE, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)
          return  ! bail out
        end if
        if (btest(verbosity,8)) then
          call ESMF_LogWrite(trim(name)//": "//trim(pName)//": Written metadata: "&
            //trim(this % meta(item) % key), ESMF_LOGMSG_INFO, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) &
            return  ! bail out
        end if
      end do

      ! close file
      call this % io % close()
      if (this % io % err % check(msg="Failure closing file "//trim(fileName), &
        line=__LINE__,  &
        file=__FILE__)) then
        call ESMF_LogSetError(ESMF_RC_FILE_CLOSE, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__,  &
          rcToReturn=rc)
        return  ! bail out
      end if
      if (btest(verbosity,8)) then
        call ESMF_LogWrite(trim(name)//": "//trim(pName)//": Closed file "//trim(fileName), &
          ESMF_LOGMSG_INFO, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
      end if

    else
      if (btest(verbosity,8)) then
        call ESMF_LogWrite(trim(name)//": "//trim(pName)//": Field: "&
          //trim(standardNameList(item))//" - time is invalid", &
          ESMF_LOGMSG_INFO, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
      end if
    end if

    if (associated(standardNameList)) then
      deallocate(standardNameList, stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
        msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
    end if

  end subroutine FileWrite

  subroutine MetadataParse(gcomp, label, phaseName, rc)
    type(ESMF_GridComp)                     :: gcomp
    character(len=*),           intent(in)  :: label
    character(len=*), optional, intent(in)  :: phaseName
    integer,          optional, intent(out) :: rc

    ! local variables
    logical                             :: isPresent
    logical                             :: eolFlag, listEnd
    integer                             :: localrc, stat
    integer                             :: item
    integer                             :: lineCount, columnCount
    integer                             :: verbosity
    character(len=ESMF_MAXSTR)          :: name
    character(len=ESMF_MAXSTR)          :: pName
    character(len=ESMF_MAXSTR)          :: msgString
    type(ESMF_Config)                   :: config
    type(SWIO_InternalState_T)          :: is
    type(SWIO_Data_T), pointer          :: this

    ! begin
    if (present(rc)) rc = ESMF_SUCCESS

    pName = "MetadataParse"
    if (present(phaseName)) pName = phaseName

    ! get component's info
    call NUOPC_CompGet(gcomp, name=name, verbosity=verbosity, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    ! get component's internal state
    call ESMF_GridCompGetInternalState(gcomp, is, localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out
    this => is % wrap

    if (this % fieldCount == 0) then
      if (btest(verbosity,8)) then
        call ESMF_LogWrite(trim(name)//": "//trim(pName)//": No input fields", &
          ESMF_LOGMSG_INFO, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
      end if
      return
    end if

    ! get component's configuration
    call ESMF_GridCompGet(gcomp, config=config, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    ! look for compute field table
    call ESMF_ConfigFindLabel(config, label, isPresent=isPresent, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    item = 0

    if (isPresent) then

      ! - get number of rows in the table
      call ESMF_ConfigGetDim(config, lineCount, columnCount, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

      ! - create metadata table in memory
      allocate(this % meta(lineCount), stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, &
        msg="Unable to allocate memory", &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

      ! - reposition input pointer to beginning of table
      call ESMF_ConfigFindLabel(config, label, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

      do
        ! get next row
        call ESMF_ConfigNextLine(config, tableEnd=listEnd, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
        if (listEnd) exit

        item = item + 1

        this % meta(item) % key   = ""
        this % meta(item) % value = ""

        ! get key
        call ESMF_ConfigGetAttribute(config, this % meta(item) % key, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out

        ! get value
        call ESMF_ConfigGetAttribute(config, this % meta(item) % value, &
          default="N/A", eolFlag=eolFlag, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out

        if (btest(verbosity,8)) then
          write(msgString,'(a,": ",a,": metadata[",i0,"]: ",a,": ",a)') &
            trim(name), trim(pName), item, trim(this % meta(item) % key), &
            trim(this % meta(item) % value)
          call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) &
            return  ! bail out
        end if
      end do

    end if

    if (btest(verbosity,8)) then
      if (item == 0) then
        call ESMF_LogWrite(trim(name)//": "//trim(pName)//": metadata: None", &
          ESMF_LOGMSG_INFO, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out
      end if
    end if

  end subroutine MetadataParse


  subroutine OutputFieldsParse(gcomp, label, phaseName, rc)
    type(ESMF_GridComp)                     :: gcomp
    character(len=*),           intent(in)  :: label
    character(len=*), optional, intent(in)  :: phaseName
    integer,          optional, intent(out) :: rc

    ! local variables
    logical                             :: isPresent
    logical                             :: eolFlag, listEnd
    integer                             :: localrc, stat
    integer                             :: item
    integer                             :: lineCount, columnCount
    integer                             :: verbosity
    character(len=ESMF_MAXSTR)          :: name
    character(len=ESMF_MAXSTR)          :: pName
    character(len=ESMF_MAXSTR)          :: msgString
    type(ESMF_Config)                   :: config
    type(SWIO_InternalState_T)          :: is
    type(SWIO_Data_T), pointer          :: this

    ! begin
    if (present(rc)) rc = ESMF_SUCCESS

    pName = "OutputFieldsParse"
    if (present(phaseName)) pName = phaseName

    ! get component's info
    call NUOPC_CompGet(gcomp, name=name, verbosity=verbosity, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    ! get component's internal state
    call ESMF_GridCompGetInternalState(gcomp, is, localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out
    this => is % wrap

    if (this % fieldCount == 0) then
      if (btest(verbosity,8)) then
        call ESMF_LogWrite(trim(name)//": "//trim(pName)//": No input fields", &
          ESMF_LOGMSG_INFO, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
      end if
      return
    end if

    ! get component's configuration
    call ESMF_GridCompGet(gcomp, config=config, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    ! look for compute field table
    call ESMF_ConfigFindLabel(config, label, isPresent=isPresent, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    item = 0

    if (isPresent) then

      ! - get number of rows in the table
      call ESMF_ConfigGetDim(config, lineCount, columnCount, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

      ! - create metadata table in memory
      allocate(this % output(lineCount), stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, &
        msg="Unable to allocate memory", &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

      ! - reposition input pointer to beginning of table
      call ESMF_ConfigFindLabel(config, label, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

      do
        ! get next row
        call ESMF_ConfigNextLine(config, tableEnd=listEnd, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
        if (listEnd) exit

        item = item + 1

        this % output(item) % key   = ""
        this % output(item) % value = ""

        ! get key
        call ESMF_ConfigGetAttribute(config, this % output(item) % key, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out

        ! get value
        call ESMF_ConfigGetAttribute(config, this % output(item) % value, &
          default=this % output(item) % key, eolFlag=eolFlag, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out

        if (btest(verbosity,8)) then
          write(msgString,'(a,": ",a,": output[",i0,"]: ",a,": ",a)') &
            trim(name), trim(pName), item, trim(this % output(item) % key), &
            trim(this % output(item) % value)
          call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) &
            return  ! bail out
        end if
      end do

    end if

    this % outputCount = this % outputCount + item

    if (btest(verbosity,8)) then
      if (item == 0) then
        call ESMF_LogWrite(trim(name)//": "//trim(pName)//": output: None", &
          ESMF_LOGMSG_INFO, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out
      end if
    end if

  end subroutine OutputFieldsParse

  subroutine OutputMaskParse(gcomp, label, phaseName, rc)
    type(ESMF_GridComp)                     :: gcomp
    character(len=*),           intent(in)  :: label
    character(len=*), optional, intent(in)  :: phaseName
    integer,          optional, intent(out) :: rc

    ! local variables
    logical                             :: isPresent
    integer                             :: localrc, stat
    integer                             :: verbosity
    integer                             :: dimCount
    character(len=ESMF_MAXSTR)          :: name
    character(len=ESMF_MAXSTR)          :: pName
    character(len=ESMF_MAXSTR)          :: fieldName
    character(len=ESMF_MAXSTR)          :: msgString
    type(ESMF_Config)                   :: config
    type(ESMF_Field)                    :: field
    type(ESMF_Grid)                     :: grid
    type(ESMF_State)                    :: importState
    type(SWIO_InternalState_T)          :: is
    type(SWIO_Data_T), pointer          :: this

    ! begin
    if (present(rc)) rc = ESMF_SUCCESS

    pName = "OutputMaskParse"
    if (present(phaseName)) pName = phaseName

    ! get component's info
    call NUOPC_CompGet(gcomp, name=name, verbosity=verbosity, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    ! get component's internal state
    call ESMF_GridCompGetInternalState(gcomp, is, localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out
    this => is % wrap

    if (trim(this % gridType) == "none") then
      if (btest(verbosity,8)) then
        call ESMF_LogWrite(trim(name)//": "//trim(pName)//&
          ": Masking only available for local 2D grids", &
          ESMF_LOGMSG_WARNING, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
      end if
      return
    end if

    if (this % fieldCount == 0) then
      if (btest(verbosity,8)) then
        call ESMF_LogWrite(trim(name)//": "//trim(pName)//": No input fields", &
          ESMF_LOGMSG_INFO, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
      end if
      return
    end if

    ! get component's configuration
    call ESMF_GridCompGet(gcomp, config=config, importState=importState, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    ! look for compute field table
    call ESMF_ConfigFindLabel(config, label, isPresent=isPresent, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    if (isPresent) then

      ! - get mask field name
      call ESMF_ConfigGetAttribute(config, fieldName, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

      ! -- get field
      call ESMF_StateGet(importState, trim(fieldName), field, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

      call ESMF_FieldGet(field, grid=grid, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

      call ESMF_GridGet(grid, dimCount=dimCount, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

      if (dimCount /= 2) return

      allocate(this % mask, stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, &
        msg="Unable to allocate memory", &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

      this % mask % field = field

      ! - get mask value
      call ESMF_ConfigGetAttribute(config, this % mask % value, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

      ! - get fill value
      call ESMF_ConfigGetAttribute(config, this % mask % fill, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

      if (btest(verbosity,8)) then
        write(msgString,'(a,": ",a,": masking: fill with ",g0.1," where ",a," == ",g0.1)') &
          trim(name), trim(pName), this % mask % fill, trim(fieldName), this % mask % value
        call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
      end if

    end if

  end subroutine OutputMaskParse

end module swio_methods
