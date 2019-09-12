c
C#############################################################################################
c
      SUBROUTINE GetNFfromVariable(Variable,NF,Variable1)
c
C#############################################################################################
c
      use User0, only: NumberOfScalarsToSolve,ScalarName,rFieldName,
     *                 NumberOfrFieldsToSolve
      use Scalar2, only: iScalarVariable
      use VolumeOfFluid2, only: irFieldVariable
c********************************************************************************************
      implicit none
c********************************************************************************************
      character*10 Variable
      character*35 Variable1
      integer :: NF,NF1
c********************************************************************************************
c
      if(Variable.eq.'velx') then
c
        NF=1
        Variable1='uvelocity'
c
      elseif(Variable.eq.'vely') then
c
        NF=2
        Variable1='vvelocity'
c
      elseif(Variable.eq.'velz') then
c
        NF=3
        Variable1='wvelocity'
c
      elseif(Variable.eq.'pressc') then
c
        NF=4
        Variable1='pressurecorrection'
c
      elseif(Variable.eq.'press') then
c
        NF=5
        Variable1='pressure'
c
      elseif(Variable.eq.'temp') then
c
        NF=6
        Variable1='temperature'
c
      elseif(Variable.eq.'htotal') then
c
        NF=7
        Variable1='temperature'
c
      elseif(Variable.eq.'tke') then
c
        NF=8
        Variable1='turbulencekineticenergy'
c
      elseif(Variable.eq.'ted') then
c
        NF=9
        Variable1='turbulencedissipationrate'
c
      elseif(Variable.eq.'tomega') then
c
        NF=10
        Variable1='turbulencespecificdissipationrate'
c
      elseif(Variable.eq.'med') then
c
        NF=11
        Variable1='modifiededdydiffusivity'
c
      elseif(Variable.eq.'tkl') then
c
        NF=12
        Variable1='turbulentkl'
c
      elseif(Variable.eq.'tgamma') then
c
        NF=13
        Variable1='turbulentgamma'
c
      elseif(Variable.eq.'tretheta') then
c
        NF=14
        Variable1='turbulentretheta'
c
      elseif(Variable.eq.'tv2') then
c
        NF=15
        Variable1='turbulentv2'
c
      elseif(Variable.eq.'tzeta') then
c
        NF=16
        Variable1='turbulentzeta'
c
      elseif(Variable.eq.'frelax') then
c
        NF=17
        Variable1='tfrelaxation'
c
      elseif(Variable.eq.'lambda') then
c
        NF=18
        Variable1='lambdaele'
c
      endif
c      
      if(irFieldVariable.ne.0) then
        if(Variable.eq.rFieldName(irFieldVariable)) then
c
          NF=18+irFieldVariable
          Variable1=rFieldName(irFieldVariable)
c
        endif
      endif
c
      if(iScalarVariable.ne.0) then
        NF1=18+NumberOfrFieldsToSolve
c      
        if(Variable.eq.ScalarName(iScalarVariable)) then
c
          NF=NF1+iScalarVariable
          Variable1=ScalarName(iScalarVariable)
c
        endif
      endif
c
      return
      end
