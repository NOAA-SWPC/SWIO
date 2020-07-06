module swio_calculator

  use ESMF
  use NUOPC
  use swio_data

  implicit none

  type SWIO_Math_T
    character(ESMF_MAXSTR) :: funcName 
    integer                :: argCount
    integer                :: parCount
    integer                :: resCount
    integer                :: rank
  end type

  type(SWIO_Math_T), parameter :: mathTable(5) = &
    (/ &
      SWIO_Math_T( "column_integrate",   1, 0, 1, 2 ),  &
      SWIO_Math_T( "column_max_point",   1, 0, 2, 2 ),  &
      SWIO_Math_T( "column_max_region",  1, 2, 2, 2 ),  &
      SWIO_Math_T( "column_interpolate", 2, 1, 1, 2 ),  &
      SWIO_Math_T( "column_o_n2_ratio",  4, 1, 1, 2 )   &
    /)

  private

  public :: SWIO_CalculatorRun
  public :: SWIO_CalculatorParse

  interface SWIO_CalculatorParse
    module procedure CalculatorParse
  end interface

  interface SWIO_CalculatorRun
    module procedure CalculatorRun
  end interface


contains

  subroutine CalculatorParse(gcomp, label, phaseName, rc)
    type(ESMF_GridComp)                     :: gcomp
    character(len=*),           intent(in)  :: label
    character(len=*), optional, intent(in)  :: phaseName
    integer,          optional, intent(out) :: rc

    ! local variables
    logical                             :: isPresent
    logical                             :: listEnd
    integer                             :: localrc, stat
    integer                             :: i, iop, item, j
    integer                             :: lineCount, columnCount
    integer                             :: verbosity
    character(len=ESMF_MAXSTR)          :: name
    character(len=ESMF_MAXSTR)          :: pName
    character(len=ESMF_MAXSTR)          :: fieldName, fieldUnits
    character(len=ESMF_MAXSTR)          :: svalue
    character(len=ESMF_MAXSTR)          :: msgString
    character(len=ESMF_MAXSTR), pointer :: standardNameList(:)
    type(ESMF_Config)                   :: config
    type(ESMF_Field)                    :: field
    type(ESMF_State)                    :: importState
    type(SWIO_InternalState_T)          :: is
    type(SWIO_Data_T), pointer          :: this

    ! begin
    if (present(rc)) rc = ESMF_SUCCESS

    pName = "CalculatorParse"
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
    call ESMF_GridCompGet(gcomp, config=config, importState=importState, &
      rc=localrc)
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

      allocate(this % task(lineCount), stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, &
        msg="Unable to allocate memory", &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

      ! - retrieve imported field list
      nullify(standardNameList)
      call NUOPC_GetStateMemberLists(importState, &
        StandardNameList=standardNameList, nestedFlag=.true., rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

      ! - read field standard name from 1st column and advertise field

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

        this % task(item) % operation = ""
        nullify(this % task(item) % fieldInp)
        nullify(this % task(item) % paramInp)

        ! get function
        call ESMF_ConfigGetAttribute(config, svalue, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out

        this % task(item) % operation = ESMF_UtilStringLowerCase(svalue, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out

        iop = 0
        do i = 1, size(mathTable)
          if (trim(this % task(item) % operation) == trim(mathTable(i) % funcName)) then
            iop = i
            exit
          end if
        end do

        if (iop == 0) then
          call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
            msg="- operation: "//trim(this % task(item) % operation), &
            line=__LINE__, &
            file=__FILE__,  &
            rcToReturn=rc)
          return  ! bail out
        end if

        allocate(this % task(item) % fieldInp(mathTable(iop) % argCount), &
          this % task(item) % fieldOut(mathTable(iop) % resCount), stat=stat)
        if (ESMF_LogFoundAllocError(statusToCheck=stat, &
          msg="Unable to allocate memory", &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out

        ! retrieve operands
        do i = 1, size(this % task(item) % fieldInp)
          ! get input field name
          call ESMF_ConfigGetAttribute(config, svalue, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) &
            return  ! bail out
          ! retrieve field from import state and store reference to operand array
          do j = 1, size(standardNameList)
            if (trim(standardNameList(j)) == trim(svalue)) then
              call ESMF_StateGet(importState, trim(standardNameList(j)), &
                this % task(item) % fieldInp(i), rc=localrc)
              if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__,  &
                file=__FILE__,  &
                rcToReturn=rc)) &
                return  ! bail out
              exit
            end if
          end do
        end do

        ! retrieve additional parameters, if required
        if (mathTable(iop) % parCount > 0) then
          allocate(this % task(item) % paramInp(mathTable(iop) % parCount), stat=stat)
          if (ESMF_LogFoundAllocError(statusToCheck=stat, &
            msg="Unable to allocate memory", &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) &
            return  ! bail out
        end if

        do i = 1, mathTable(iop) % parCount
          call ESMF_ConfigGetAttribute(config, this % task(item) % paramInp(i), &
            rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) &
            return  ! bail out
        end do

        ! create output fields
        do i = 1, size(this % task(item) % fieldOut)
          ! get computed field name
          call ESMF_ConfigGetAttribute(config, fieldName, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) &
            return  ! bail out
          ! get computed field units
          call ESMF_ConfigGetAttribute(config, fieldUnits, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) &
            return  ! bail out
          ! create computed field
          this % task(item) % fieldOut(i) = &
            FieldCreate(this % task(item) % fieldInp(1), mathTable(iop) % rank, &
              fieldName, fieldUnits, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) &
            return  ! bail out
        end do

        ! retrieve scaling factor, if present
        call ESMF_ConfigGetAttribute(config, this % task(item) % scaleFactor, &
          default=1._ESMF_KIND_R8, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out

        if (btest(verbosity,8)) then
          write(msgString,'(a,": ",a,": compute[",i0,"]: ")') trim(name), &
            trim(pName), item
          call AppendArgumentListString(msgString, this % task(item) % fieldOut, &
            rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) &
            return  ! bail out
          call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) &
            return  ! bail out
          write(msgString,'(a,": ",a,": - Function  : ",a)') trim(name), &
            trim(pName), trim(this % task(item) % operation)
          call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) &
            return  ! bail out
          write(msgString,'(a,": ",a,": - Arguments : ")') trim(name), trim(pName)
          call AppendArgumentListString(msgString, this % task(item) % fieldInp, &
            paramList=this % task(item) % paramInp, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) &
            return  ! bail out
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
        call ESMF_LogWrite(trim(name)//": "//trim(pName)//": compute: None", &
          ESMF_LOGMSG_INFO, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out
      end if
    end if

  end subroutine CalculatorParse


  subroutine Calculate(task, rc)
    type(SWIO_Task_T), intent(inout) :: task
    integer, optional, intent(out)   :: rc

    ! -- local variables
    integer :: localrc

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    select case (trim(task % operation))
      case ("column_integrate")
        call columnIntegrate(task, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, &
          msg="Failure in function: "//task % operation, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
      case ("column_max_point")
        call columnMaxLoc(task, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, &
          msg="Failure in function: "//task % operation, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
      case ("column_max_region")
        call columnMaxLocRegion(task, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, &
          msg="Failure in function: "//task % operation, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
      case ("column_interpolate")
        call columnInterpolate(task, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, &
          msg="Failure in function: "//task % operation, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
      case ("column_o_n2_ratio")
        call columnON2ratio(task, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, &
          msg="Failure in function: "//task % operation, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
      case default
        ! -- nothing to do
    end select

  end subroutine Calculate

  subroutine CalculatorRun(gcomp, phaseName, rc)
    type(ESMF_GridComp)                     :: gcomp
    character(len=*), optional, intent(in)  :: phaseName
    integer,          optional, intent(out) :: rc

    ! local variables
    logical                             :: isPresent
    integer                             :: localrc
    integer                             :: i, iop, item, stat
    integer                             :: verbosity
    character(len=ESMF_MAXSTR)          :: name
    character(len=ESMF_MAXSTR)          :: msgString
    character(len=ESMF_MAXSTR)          :: pName
    character(len=ESMF_MAXSTR), pointer :: standardNameList(:)
    type(ESMF_Field)                    :: field
    type(SWIO_InternalState_T)          :: is
    type(SWIO_Data_T), pointer          :: this

    ! begin
    if (present(rc)) rc = ESMF_SUCCESS

    pName = "CalculatorRun"
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

    ! run calculator's tasks
    do item = 1, size(this % task)
      call Calculate(this % task(item), rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      if (btest(verbosity,8)) then
          write(msgString,'(a,": computed[",i0,"]: ")') trim(name) &
            //": "//trim(pName), item
          call AppendArgumentListString(msgString, this % task(item) % fieldOut, &
            rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) &
            return  ! bail out
          call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) &
            return  ! bail out
      end if
    end do

  end subroutine CalculatorRun

  ! -- Tools
  subroutine AppendArgumentListString(string, fieldList, paramList, rc)
    character(len=*),             intent(inout) :: string
    type(ESMF_Field),             pointer       :: fieldList(:)
    real(ESMF_KIND_R8), optional, pointer       :: paramList(:)
    integer,            optional, intent(out)   :: rc

    ! -- local variables
    integer :: localrc
    integer :: i, item
    character(ESMF_MAXSTR) :: fieldName

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    item = 0
    if (associated(fieldList)) then
      do i = 1, size(fieldList)
        call ESMF_FieldGet(fieldList(i), name=fieldName, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
        item = item + 1
        if (item > 1) then
          string = trim(string) // ", " // trim(fieldName)
        else
          string = trim(string) // "  " // trim(fieldName)
        end if
      end do
    end if

    if (present(paramList)) then
      if (associated(paramList)) then
        do i = 1, size(paramList)
          write(fieldName, '(g0.5)', iostat=localrc) paramList(i)
          if (localrc /= 0) then
            call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__, &
              file=__FILE__, &
              rcToReturn=rc)
            return  ! bail out
          end if
          item = item + 1
          if (item > 1) then
            string = trim(string) // ", " // trim(fieldName)
          else
            string = trim(string) // "  " // trim(fieldName)
          end if
        end do
      end if
    end if

  end subroutine AppendArgumentListString

  function FieldCreate(field, rank, name, units, rc) result (fieldOut)
    type(ESMF_Field),         intent(in)  :: field
    integer,                  intent(in)  :: rank
    character(len=*),         intent(in)  :: name
    character(len=*),         intent(in)  :: units
    integer, optional,        intent(out) :: rc

    integer                  :: localrc, stat
    integer                  :: dimCount
    integer, allocatable     :: gridToFieldMap(:)
    character(ESMF_MAXSTR)   :: fieldName
    type(ESMF_Field)         :: fieldOut
    type(ESMF_Grid)          :: grid
    type(ESMF_GeomType_Flag) :: geomType
    type(ESMF_Mesh)          :: mesh
    type(ESMF_MeshLoc)       :: meshloc
    type(ESMF_StaggerLoc)    :: staggerloc
    type(ESMF_TypeKind_Flag) :: typekind

    ! begin
    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldGet(field, name=fieldName, geomtype=geomType, &
      typekind=typekind, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    if (geomType == ESMF_GEOMTYPE_GRID) then

      call ESMF_FieldGet(field, grid=grid, staggerloc=staggerloc, rc=localrc)
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

      if (rank > dimCount) then
        ! -- add ungridded dimensions (not yet implemented)
        call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
          msg="output field rank must be <= input field dimCount", &
          line=__LINE__, &
          file=__FILE__,  &
          rcToReturn=rc)
        return  ! bail out
      else
        allocate(gridToFieldMap(dimCount), stat=stat)
        if (ESMF_LogFoundAllocError(statusToCheck=stat, &
          msg="Unable to allocate memory", &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
        call ESMF_FieldGet(field, gridToFieldMap=gridToFieldMap, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
        if (rank < dimCount) gridToFieldMap(rank+1:) = 0 
      end if

      fieldOut = ESMF_FieldCreate(grid, typekind, staggerloc=staggerloc, &
        gridToFieldMap=gridToFieldMap, name=name, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

      deallocate(gridToFieldMap, stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
        msg="Unable to free up memory", &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

    else if (geomType == ESMF_GEOMTYPE_MESH) then

      call ESMF_FieldGet(field, mesh=mesh, meshloc=meshloc, rc=localrc)
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

      if (rank > dimCount) then
        ! -- add ungridded dimensions (not yet implemented)
        call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
          msg="output field rank must be <= input field dimCount", &
          line=__LINE__, &
          file=__FILE__,  &
          rcToReturn=rc)
        return  ! bail out
      else
        allocate(gridToFieldMap(dimCount), stat=stat)
        if (ESMF_LogFoundAllocError(statusToCheck=stat, &
          msg="Unable to allocate memory", &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
        call ESMF_FieldGet(field, gridToFieldMap=gridToFieldMap, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) &
          return  ! bail out
        if (rank < dimCount) gridToFieldMap(rank+1:) = 0 
      end if

      fieldOut = ESMF_FieldCreate(mesh, typekind, meshloc=meshloc, &
        gridToFieldMap=gridToFieldMap, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

      deallocate(gridToFieldMap, stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
        msg="Unable to free up memory", &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

    else
      call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
        msg="Only Grid and Mesh geometry types are supported", &
        line=__LINE__, &
        file=__FILE__,  &
        rcToReturn=rc)
      return  ! bail out
    end if

    ! add field metadata and initialize

    ! - add Attribute packages
    call ESMF_AttributeAdd(fieldOut, convention="ESG", purpose="General", &
      rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out
    call ESMF_AttributeAdd(fieldOut, convention="NUOPC", purpose="Instance",   &
      nestConvention="ESG", nestPurpose="General", rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    ! - add name and units
    call ESMF_AttributeSet(fieldOut, &
      name="StandardName", value=name, &
      convention="NUOPC", purpose="Instance", attnestflag=ESMF_ATTNEST_ON, &
      rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out
    call ESMF_AttributeSet(fieldOut, &
      name="LongName", value=name, &
      convention="NUOPC", purpose="Instance", attnestflag=ESMF_ATTNEST_ON, &
      rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out
    call ESMF_AttributeSet(fieldOut, &
      name="Units", value=units, &
      convention="NUOPC", purpose="Instance", attnestflag=ESMF_ATTNEST_ON, &
      rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    ! - initialize to zero
    call ESMF_FieldFill(fieldOut, dataFillScheme="const", &
      const1=0._ESMF_KIND_R8, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

  end function FieldCreate

  logical function IsTaskValid(task, table, rc)
    type(SWIO_Task_T), intent(in)  :: task
    type(SWIO_Math_T), intent(in)  :: table
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: item

    ! -- begin
    IsTaskValid = .false.
    if (present(rc)) rc = ESMF_SUCCESS

    if (associated(task % fieldInp)) then
      if (size(task % fieldInp) /= table % argCount) then
        call ESMF_LogSetError(ESMF_RC_INTNRL_INCONS, &
          msg="wrong number of input field arguments", &
          line=__LINE__, &
          file=__FILE__,  &
          rcToReturn=rc)
        return  ! bail out
      end if
    else
      call ESMF_LogSetError(ESMF_RC_PTR_NOTALLOC, &
        msg="missing input field arguments", &
        line=__LINE__, &
        file=__FILE__,  &
        rcToReturn=rc)
      return  ! bail out
    end if

    if (associated(task % fieldOut)) then
      if (size(task % fieldOut) /= table % resCount) then
        call ESMF_LogSetError(ESMF_RC_INTNRL_INCONS, &
          msg="wrong number of output field arguments", &
          line=__LINE__, &
          file=__FILE__,  &
          rcToReturn=rc)
        return  ! bail out
      end if
    else
      call ESMF_LogSetError(ESMF_RC_PTR_NOTALLOC, &
        msg="missing output field arguments", &
        line=__LINE__, &
        file=__FILE__,  &
        rcToReturn=rc)
      return  ! bail out
    end if

    if (table % parCount > 0) then
      if (associated(task % paramInp)) then
        if (size(task % paramInp) /= table % parCount) then
          call ESMF_LogSetError(ESMF_RC_INTNRL_INCONS, &
            msg="wrong number of input parameters", &
            line=__LINE__, &
            file=__FILE__,  &
            rcToReturn=rc)
          return  ! bail out
        end if
      else
        call ESMF_LogSetError(ESMF_RC_PTR_NOTALLOC, &
          msg="missing input parameters", &
          line=__LINE__, &
          file=__FILE__,  &
          rcToReturn=rc)
        return  ! bail out
      end if
    end if

    IsTaskValid = .true.

  end function IsTaskValid


  ! -- Available math functions

  subroutine columnIntegrate(task, rc)
    type(SWIO_Task_T), intent(inout) :: task
    integer, optional, intent(out)   :: rc

    ! -- local variables
    integer :: localrc
    integer :: dimCount, localDe, localDeCount
    integer :: i, j, k, km1
    integer, dimension(3) :: lb, ub
    logical :: isValid
    real(ESMF_KIND_R8) :: coef
    real(ESMF_KIND_R8), dimension(:),     pointer :: z
    real(ESMF_KIND_R8), dimension(:,:),   pointer :: q
    real(ESMF_KIND_R8), dimension(:,:,:), pointer :: f
    type(ESMF_Grid)  :: grid

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    ! -- check if arguments are valid
    isValid = isTaskValid(task, mathTable(1), rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out
    if (.not.isValid) then
      call ESMF_LogSetError(ESMF_RC_ARG_INCOMP, &
        msg="invalid argument list", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if

    call ESMF_FieldGet(task % fieldInp(1), grid=grid, &
      localDeCount=localDeCount, rc=localrc)
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

    if (dimCount /= 3) then
      call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
        msg="this function only supports fields on 3D grids", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if

    do localDe = 0, localDeCount-1
      call ESMF_FieldGet(task % fieldInp(1), localDe=localDe, farrayPtr=f, &
        computationalLBound=lb, computationalUBound=ub, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      call ESMF_GridGetCoord(grid, 3, localDe=localDe, farrayPtr=z, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      call ESMF_FieldGet(task % fieldOut(1), localDe=localDe, farrayPtr=q, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

      q(lb(1):ub(1),lb(2):ub(2)) = 0._ESMF_KIND_R8

      do k = lb(3) + 1, ub(3)
        km1 = k - 1
        coef = 0.5_ESMF_KIND_R8 * task % scaleFactor * (z(k) - z(km1))
        do j = lb(2), ub(2)
          do i = lb(1), ub(1)
            q(i,j) = q(i,j) + coef * (f(i,j,km1)+f(i,j,k))
          end do
        end do
      end do
    end do

  end subroutine columnIntegrate

  subroutine columnMaxLoc(task, rc)
    type(SWIO_Task_T), intent(inout) :: task
    integer, optional, intent(out)   :: rc

    ! -- local variables
    integer :: localrc
    integer :: dimCount, localDe, localDeCount
    integer :: i, j, k1, k2, k3, k_max
    integer, dimension(3) :: lb, ub
    logical :: isValid
    real(ESMF_KIND_R8) :: x1, x2, x3, y1, y2, y3, a, b, c
    real(ESMF_KIND_R8), dimension(:),     pointer :: z
    real(ESMF_KIND_R8), dimension(:,:),   pointer :: q, h
    real(ESMF_KIND_R8), dimension(:,:,:), pointer :: f
    type(ESMF_Grid)  :: grid

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    ! -- check if arguments are valid
    isValid = isTaskValid(task, mathTable(2), rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out
    if (.not.isValid) then
      call ESMF_LogSetError(ESMF_RC_ARG_INCOMP, &
        msg="invalid argument list", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if

    call ESMF_FieldGet(task % fieldOut(1), grid=grid, &
      localDeCount=localDeCount, rc=localrc)
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

    if (dimCount /= 3) then
      call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
        msg="this function only supports fields on 3D grids", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if

    do localDe = 0, localDeCount-1
      call ESMF_FieldGet(task % fieldInp(1), localDe=localDe, farrayPtr=f, &
        computationalLBound=lb, computationalUBound=ub, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      call ESMF_GridGetCoord(grid, 3, localDe=localDe, farrayPtr=z, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      call ESMF_FieldGet(task % fieldOut(1), localDe=localDe, farrayPtr=q, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      call ESMF_FieldGet(task % fieldOut(2), localDe=localDe, farrayPtr=h, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

      q(lb(1):ub(1),lb(2):ub(2)) = 0._ESMF_KIND_R8
      h(lb(1):ub(1),lb(2):ub(2)) = 0._ESMF_KIND_R8

      do j = lb(2), ub(2)
        do i = lb(1), ub(1)
          k_max = maxval(maxloc(f(i,j,:)))
          k1 = max(k_max - 1, lb(3))
          k2 = k_max
          k3 = min(k_max + 1, ub(3))
          x1 = z(k1)
          x2 = z(k2)
          x3 = z(k3)
          y1 = f(i,j,k1)
          y2 = f(i,j,k2)
          y3 = f(i,j,k3)
          c = (x3*y1 - x3*y2 + x1*y2 + x2*y3 - x2*y1 - x1*y3) &
              / (x3*x1*x1 - x3*x2*x2 + x1*x2*x2 - x2*x1*x1 + x2*x3*x3 - x1*x3*x3)
          b = (y2 - (c*x2*x2) + (c*x1*x1) - y1) / (x2 - x1)
          a = y1 - (b*x1) - (c*x1*x1)
          h(i,j) = (0._ESMF_KIND_R8 - b) / (2*c)
          q(i,j) = a + (b*h(i,j)) + (c*h(i,j)*h(i,j))
        end do
      end do

    end do

  end subroutine columnMaxLoc

  subroutine columnMaxLocRegion(task, rc)
    type(SWIO_Task_T), intent(inout) :: task
    integer, optional, intent(out)   :: rc

    ! -- local variables
    integer :: localrc
    integer :: dimCount, localDe, localDeCount
    integer :: i, j
    integer :: ka, kb, km, km1, kp1, kloc, kglb
    integer, dimension(3) :: lb, ub
    logical :: isValid
    real(ESMF_KIND_R8) :: x1, x2, x3, y1, y2, y3, a, b, c, fm, fmax
    real(ESMF_KIND_R8), dimension(:),     pointer :: z
    real(ESMF_KIND_R8), dimension(:,:),   pointer :: q, h
    real(ESMF_KIND_R8), dimension(:,:,:), pointer :: f
    type(ESMF_Grid)  :: grid

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    ! -- check if arguments are valid
    isValid = isTaskValid(task, mathTable(3), rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out
    if (.not.isValid) then
      call ESMF_LogSetError(ESMF_RC_ARG_INCOMP, &
        msg="invalid argument list", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if

    call ESMF_FieldGet(task % fieldOut(1), grid=grid, &
      localDeCount=localDeCount, rc=localrc)
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

    if (dimCount /= 3) then
      call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
        msg="this function only supports fields on 3D grids", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if

    do localDe = 0, localDeCount-1
      call ESMF_FieldGet(task % fieldInp(1), localDe=localDe, farrayPtr=f, &
        computationalLBound=lb, computationalUBound=ub, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      call ESMF_GridGetCoord(grid, 3, localDe=localDe, farrayPtr=z, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      call ESMF_FieldGet(task % fieldOut(1), localDe=localDe, farrayPtr=q, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      call ESMF_FieldGet(task % fieldOut(2), localDe=localDe, farrayPtr=h, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

      q(lb(1):ub(1),lb(2):ub(2)) = 0._ESMF_KIND_R8
      h(lb(1):ub(1),lb(2):ub(2)) = 0._ESMF_KIND_R8

      do j = lb(2), ub(2)
        do i = lb(1), ub(1)
          ka = minloc(z, dim=1, mask = z >= task % paramInp(1))
          kb = maxloc(z, dim=1, mask = z <= task % paramInp(2))

          kloc = 0
          kglb = kb
          km1  = 0
          kp1  = 0
          fmax = -1._ESMF_KIND_R8
          do km = kb-1, ka+1, -1
            km1 = max(km-1,ka)
            kp1 = min(km+1,kb)
            fm  = f(i,j,km)
            if (fm > f(i,j,kp1)) then
              kglb = km
              if (fm > f(i,j,km1) .and. fm > fmax) then
                fmax = fm
                kloc = km
              end if
            end if
          end do

          if (kloc < 1) then
            ! -- local maximum not found: use global maximum
            if (f(i,j,ka) > f(i,j,kglb)) kglb = ka
            h(i,j) = z(kglb)
            q(i,j) = f(i,j,kglb)
          else
            km1 = max(kloc-1,ka)
            kp1 = min(kloc+1,kb)
            if (kp1-km1 == 2) then
              ! -- refine maximum through quadratic fit
              x1 = z(km1)
              x2 = z(kloc)
              x3 = z(kp1)
              y1 = f(i,j,km1)
              y2 = f(i,j,kloc)
              y3 = f(i,j,kp1)
              c = (x3*y1 - x3*y2 + x1*y2 + x2*y3 - x2*y1 - x1*y3) &
                  / (x3*x1*x1 - x3*x2*x2 + x1*x2*x2 - x2*x1*x1 + x2*x3*x3 - x1*x3*x3)
              b = (y2 - (c*x2*x2) + (c*x1*x1) - y1) / (x2 - x1)
              a = y1 - (b*x1) - (c*x1*x1)
              h(i,j) = - b / (2*c)
              q(i,j) = a + (b*h(i,j)) + (c*h(i,j)*h(i,j))
            else
              ! -- not enough points for quadratic fit
              h(i,j) = z(kloc)
              q(i,j) = f(i,j,kloc)
            end if
          end if
        end do
      end do

    end do

  end subroutine columnMaxLocRegion

  subroutine columnInterpolate(task, rc)
    type(SWIO_Task_T), intent(inout) :: task
    integer, optional, intent(out)   :: rc

    ! -- local variables
    integer :: localrc
    integer :: dimCount, rank
    integer :: localDe, localDeCount
    integer :: i, j
    integer, dimension(3) :: lb, ub
    logical :: isValid
    real(ESMF_KIND_R8), dimension(1)              :: xd, yd
    real(ESMF_KIND_R8), dimension(:),     pointer :: x, y
    real(ESMF_KIND_R8), dimension(:,:),   pointer :: q
    real(ESMF_KIND_R8), dimension(:,:,:), pointer :: f, h
    type(ESMF_Grid)  :: grid

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    ! -- check if arguments are valid
    isValid = isTaskValid(task, mathTable(4), rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out
    if (.not.isValid) then
      call ESMF_LogSetError(ESMF_RC_ARG_INCOMP, &
        msg="invalid argument list", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if

    do i = 1, size(task % fieldInp)
      call ESMF_FieldGet(task % fieldInp(i), grid=grid, dimCount=dimCount, &
        rank=rank, localDeCount=localDeCount, rc=localrc)
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
      if ((dimCount /= 2) .or. (rank /= 3)) then
        call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
          msg="this function only supports 2D fields with multiple vertical levels", &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)
        return  ! bail out
      end if
    end do

    do localDe = 0, localDeCount-1
      call ESMF_FieldGet(task % fieldInp(1), localDe=localDe, farrayPtr=f, &
        computationalLBound=lb, computationalUBound=ub, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      call ESMF_FieldGet(task % fieldInp(2), localDe=localDe, farrayPtr=h, &
        rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      call ESMF_FieldGet(task % fieldOut(1), localDe=localDe, farrayPtr=q, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

      nullify(x, y)

      xd(1) = task % paramInp(1)
      q(lb(1):ub(1),lb(2):ub(2)) = 0._ESMF_KIND_R8

      do j = lb(2), ub(2)
        do i = lb(1), ub(1)
          x => h(i,j,:)
          y => f(i,j,:)
          yd(1) = 0._ESMF_KIND_R8
          call PolyInterpolate(x, y, xd, yd, 1, localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) &
            return  ! bail out
          q(i,j) = yd(1)
        end do
      end do

    end do

  end subroutine columnInterpolate

  subroutine columnON2ratio(task, rc)
    type(SWIO_Task_T), intent(inout) :: task
    integer, optional, intent(out)   :: rc

    ! -- local variables
    integer :: localrc
    integer :: dimCount, rank
    integer :: localDe, localDeCount
    integer :: i, j, k
    integer, dimension(3) :: lb, ub
    logical :: isValid
    real(ESMF_KIND_R8) :: fscale
    real(ESMF_KIND_R8) :: n2_lev_ratio, n2_target, n2_total, n2_total_prev
    real(ESMF_KIND_R8) :: o_total, o_total_prev
    real(ESMF_KIND_R8), dimension(:,:),   pointer :: r
    real(ESMF_KIND_R8), dimension(:,:,:), pointer :: o, n2, t, h
    type(ESMF_Grid)  :: grid

    ! -- constants
    ! From: CODATA Recommended Values of the Fundamental Physical Constants: 2018
    ! - standard acceleration of gravity
    real(ESMF_KIND_R8), parameter :: g0 = 9.80665_ESMF_KIND_R8
    ! - unified atomic mass unit (kg)
    real(ESMF_KIND_R8), parameter :: um = 1.66053906660E-27_ESMF_KIND_R8
    ! - Boltzmann's constant (J K-1)
    real(ESMF_KIND_R8), parameter :: kB = 1.380649E-23_ESMF_KIND_R8
    ! From: Geodetic Reference System 1980.
    ! - Earth's mean radius (m)
    real(ESMF_KIND_R8), parameter :: Re = 6371008.7714_ESMF_KIND_R8
    ! From: Atomic weights of the elements 2013 (IUPAC Technical Report)
    !       Pure and Applied Chemistry, Volume 88, Issue 3, Pages 265–291, eISSN
    !       1365-3075, ISSN 0033-4545, DOI: https://doi.org/10.1515/pac-2015-0305.
    ! - Standard atomic weight of oxygen (air)
    real(ESMF_KIND_R8), parameter :: o_weight = 15.9994_ESMF_KIND_R8
    ! - Standard molecular weight of nitrogen (air)
    real(ESMF_KIND_R8), parameter :: n2_weight = 2 * 14.0067_ESMF_KIND_R8

    ! -- local constants
    real(ESMF_KIND_R8), parameter :: const = kB / (um * g0)

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    ! -- check if arguments are valid
    isValid = isTaskValid(task, mathTable(5), rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out
    if (.not.isValid) then
      call ESMF_LogSetError(ESMF_RC_ARG_INCOMP, &
        msg="invalid argument list", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    end if

    do i = 1, size(task % fieldInp)
      call ESMF_FieldGet(task % fieldInp(i), grid=grid, dimCount=dimCount, &
        rank=rank, localDeCount=localDeCount, rc=localrc)
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
      if ((dimCount /= 2) .or. (rank /= 3)) then
        call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
          msg="this function only supports 2D fields with multiple vertical levels", &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)
        return  ! bail out
      end if
    end do

    do localDe = 0, localDeCount-1
      call ESMF_FieldGet(task % fieldInp(1), localDe=localDe, farrayPtr=o, &
        computationalLBound=lb, computationalUBound=ub, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      call ESMF_FieldGet(task % fieldInp(2), localDe=localDe, farrayPtr=n2, &
        rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      call ESMF_FieldGet(task % fieldInp(3), localDe=localDe, farrayPtr=t, &
        rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      call ESMF_FieldGet(task % fieldInp(4), localDe=localDe, farrayPtr=h, &
        rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out
      call ESMF_FieldGet(task % fieldOut(1), localDe=localDe, farrayPtr=r, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) &
        return  ! bail out

      n2_target = task % paramInp(1)
      n2_total  = 0._ESMF_KIND_R8
      o_total   = 0._ESMF_KIND_R8
      r(lb(1):ub(1),lb(2):ub(2)) = 0._ESMF_KIND_R8

      if (n2_target > 0._ESMF_KIND_R8) then
        do j = lb(2), ub(2)
          do i = lb(1), ub(1)
            do k = lb(3), ub(3)
              fscale = 1._ESMF_KIND_R8 + h(i,j,k) / Re
              fscale = const * fscale * fscale * t(i,j,k)

              o_total_prev  = o_total
              n2_total_prev = n2_total
              n2_total = n2(i,j,k) * fscale / n2_weight
              o_total  = o (i,j,k) * fscale / o_weight

              if ((n2_total < n2_target) .and. (k > lb(3)) .and. &
                  (n2_total > 0._ESMF_KIND_R8) .and. (n2_total_prev /= n2_total) .and. &
                  ( o_total > 0._ESMF_KIND_R8) .and. ( o_total_prev /= o_total)) then
                n2_lev_ratio = log(n2_target/n2_total) / log(n2_total_prev/n2_total)
                r(i,j) = o_total * exp(n2_lev_ratio*log(o_total_prev/o_total)) / n2_target
                exit
              end if

            end do
          end do
        end do
      end if

    end do

  end subroutine columnON2ratio

  ! -- Shared auxiliary math functions

  subroutine PolyInterpolate(xs, ys, xd, yd, m, rc)
    real(ESMF_KIND_R8), dimension(:), intent(in)  :: xs, ys, xd
    real(ESMF_KIND_R8), dimension(:), intent(out) :: yd
    integer, intent(in)  :: m
    integer, intent(out) :: rc

    ! -- local variables
    integer :: i, j, k, n, np
    real(ESMF_KIND_R8) :: x, y, dy

    ! -- begin
    rc = ESMF_SUCCESS

    n = m + 1
    np = size(xs)
    do i = 1, size(xd)
      x = xd(i)
      y = 0._ESMF_KIND_R8
      call locate(xs, np, x, j)
      if (j == np) then
        y = ys(np)
      else if (j > 0) then
        k = min(max(j-(n-1)/2,1), np+1-n)
        call polint(xs(k:), ys(k:), n, x, y, dy, rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg="Error in polint", &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      end if
      yd(i) = y
    end do

  end subroutine PolyInterpolate

  ! -- auxiliary numerical subroutines for interpolation/extrapolation from:
  ! -- W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery,
  ! -- Numerical Recipes in Fortran 77: The Art of Scientific Computing
  ! -- (Vol. 1 of Fortran Numerical Recipes), 2nd Ed., Cambridge Univ. Press
  ! -- 1992

  subroutine locate(xx, n, x, j)
    implicit none
    integer, intent(in) :: n
    real(ESMF_KIND_R8), intent(in) :: x, xx(n)
    integer, intent(out) :: j

    ! -- local variables
    integer :: jl, jm, ju

    ! -- begin
    jl = 0
    ju = n + 1
    do while (ju-jl.gt.1)
      jm = (ju+jl)/2
      if ((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm))) then
        jl = jm
      else
        ju = jm
      end if
    end do
    if (x.eq.xx(1)) then
      j = 1
    else if (x.eq.xx(n)) then
      j = n - 1
    else
      j = jl
    end if

  end subroutine locate

  subroutine polint(xa, ya, n, x, y, dy, rc)

    implicit none

    integer, intent(in) :: n
    real(ESMF_KIND_R8),    intent(in) :: xa(n), ya(n)
    real(ESMF_KIND_R8),    intent(in) :: x
    real(ESMF_KIND_R8),   intent(out) :: y, dy
    integer, intent(out) :: rc

    ! -- local variables
    integer, parameter :: nmax = 10
    integer :: i, m, ns
    real(ESMF_KIND_R8) :: den, dif, dift, ho, hp, w
    real(ESMF_KIND_R8), dimension(nmax) :: c, d

    ! -- begin
    rc = ESMF_SUCCESS

    ns = 1
    dif = abs(x - xa(1))
    do i = 1, n
      dift = abs(x - xa(i))
      if (dift < dif) then
        ns  = i
        dif = dift
      end if
      c(i) = ya(i)
      d(i) = ya(i)
    end do
    y = ya(ns)
    ns = ns - 1
    do m = 1, n - 1
      do i = 1, n - m
        ho = xa(i)   - x
        hp = xa(i+m) - x
        w  = c(i+1)-d(i)
        den = ho - hp
        if (den == 0.) then
          rc = ESMF_FAILURE
          exit
        end if
        den = w / den
        d(i) = hp * den
        c(i) = ho * den
      end do
      if (2*ns < n - m) then
        dy = c(ns+1)
      else
        dy = d(ns)
        ns = ns - 1
      end if
      y = y + dy
    end do

  end subroutine polint

end module swio_calculator
