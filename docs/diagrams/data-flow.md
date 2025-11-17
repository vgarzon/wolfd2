# Data Flow Diagram

This diagram illustrates how data flows through the wolfd2 program, from input files through processing to output files.

```mermaid
flowchart LR
    subgraph Input["Input Data Sources"]
        InputFile[input.dat<br/>Simulation Parameters]
        BCFile[bcond.dat<br/>Boundary Conditions]
        GridFile[Grid File<br/>Plot3D Format]
        RestartFile[Restart File<br/>Optional]
    end
    
    subgraph Parse["Parsing Module<br/>parse.f"]
        ReadParam[ReadParam<br/>Parse Parameters]
        ReadBC[ReadBC<br/>Parse BCs]
    end
    
    subgraph Grid["Grid Module<br/>grid.f"]
        ReadGrid[ReadGrid<br/>Load Coordinates]
        Metric[Metric<br/>Compute Metrics]
    end
    
    subgraph Init["Initialization<br/>bound_cond.f"]
        SetupBC[SetUpBCs<br/>Configure BCs]
        InitialCond[InitialCond<br/>Initialize Fields]
    end
    
    subgraph SolnVars["Solution Variables"]
        U[Velocity u]
        V[Velocity v]
        P[Pressure p]
        T[Temperature t]
        D[Density d]
        USS[Small-Scale u]
        VSS[Small-Scale v]
        PSS[Small-Scale p]
        TSS[Small-Scale t]
    end
    
    subgraph GridVars["Grid Variables"]
        GX[Grid x-coords]
        GY[Grid y-coords]
        Metrics[Metric Tensors<br/>ran, rbn, rgn, etc.]
        Jacobians[Jacobians<br/>dju, djv, djc, djn]
    end
    
    subgraph BCVars["Boundary Condition Arrays"]
        RegType[Region Types]
        MomBdTp[Momentum BC Types]
        TemBdTp[Thermal BC Types]
        BCVal[BC Values]
    end
    
    subgraph Momentum["Momentum Solver<br/>momentum.f"]
        AuxMom[nAuxMomentum<br/>Auxiliary Velocities]
        XMom[XMomentum<br/>X-direction]
        YMom[YMomentum<br/>Y-direction]
    end
    
    subgraph Pressure["Pressure Solver<br/>pressure.f"]
        Div[Divergence<br/>Compute div u*]
        PPE[Ppe<br/>Solve Poisson Eq]
        SOR[SOR Solvers<br/>Iterative Solution]
    end
    
    subgraph Thermal["Thermal Solver<br/>thermal.f"]
        ThermEnergy[ThermEnergy<br/>Solve Energy Eq]
        EqState[EqState<br/>Compute Density]
    end
    
    subgraph SmallScale["Small-Scale Module<br/>small_scale.f"]
        SmallScaleFunc[SmallScale<br/>ATD Computation]
    end
    
    subgraph Utility["Utility Functions<br/>utility.f"]
        Filter[Filter<br/>Shuman Filter]
        Project[Project<br/>Velocity Projection]
        Avg[Average Functions<br/>Interpolation]
    end
    
    subgraph Output["Output Files"]
        LogFile[output.log<br/>Execution Log]
        SolnFile[output.xyz<br/>Solution Fields]
        QQQFile[output.qqq<br/>Additional Fields]
        SnapFiles[Snapshot Files<br/>Time-Dependent]
        TSFiles[Time Series Files<br/>Point Histories]
        RestartOut[Restart Files<br/>Checkpoints]
        TrajFiles[Trajectory Files<br/>Particle Paths]
    end
    
    InputFile --> ReadParam
    BCFile --> ReadBC
    GridFile --> ReadGrid
    RestartFile --> InitialCond
    
    ReadParam --> SolnVars
    ReadParam --> BCVars
    ReadBC --> BCVars
    ReadGrid --> GX
    ReadGrid --> GY
    GX --> Metric
    GY --> Metric
    Metric --> Metrics
    Metric --> Jacobians
    
    Metrics --> Momentum
    Metrics --> Pressure
    Metrics --> Thermal
    Metrics --> SmallScale
    Jacobians --> Momentum
    Jacobians --> Pressure
    Jacobians --> Thermal
    
    BCVars --> Init
    Init --> SetupBC
    SetupBC --> RegType
    SetupBC --> MomBdTp
    SetupBC --> TemBdTp
    SetupBC --> BCVal
    
    InitialCond --> U
    InitialCond --> V
    InitialCond --> P
    InitialCond --> T
    
    U --> AuxMom
    V --> AuxMom
    RegType --> AuxMom
    MomBdTp --> AuxMom
    BCVal --> AuxMom
    Metrics --> AuxMom
    AuxMom --> XMom
    AuxMom --> YMom
    XMom --> U
    YMom --> V
    
    U --> Filter
    V --> Filter
    Filter --> U
    Filter --> V
    
    U --> Div
    V --> Div
    Div --> PPE
    PPE --> SOR
    SOR --> P
    
    P --> Project
    U --> Project
    V --> Project
    Project --> U
    Project --> V
    
    U --> ThermEnergy
    V --> ThermEnergy
    T --> ThermEnergy
    RegType --> ThermEnergy
    TemBdTp --> ThermEnergy
    BCVal --> ThermEnergy
    ThermEnergy --> T
    
    P --> EqState
    T --> EqState
    EqState --> D
    
    U --> SmallScaleFunc
    V --> SmallScaleFunc
    T --> SmallScaleFunc
    RegType --> SmallScaleFunc
    BCVal --> SmallScaleFunc
    SmallScaleFunc --> USS
    SmallScaleFunc --> VSS
    SmallScaleFunc --> PSS
    SmallScaleFunc --> TSS
    
    USS --> U
    VSS --> V
    TSS --> T
    
    U --> Avg
    V --> Avg
    P --> Avg
    T --> Avg
    USS --> Avg
    VSS --> Avg
    
    U --> LogFile
    V --> LogFile
    P --> LogFile
    T --> LogFile
    U --> SolnFile
    V --> SolnFile
    P --> SolnFile
    T --> SolnFile
    D --> QQQFile
    U --> SnapFiles
    V --> SnapFiles
    P --> SnapFiles
    T --> SnapFiles
    U --> TSFiles
    V --> TSFiles
    P --> TSFiles
    T --> TSFiles
    U --> RestartOut
    V --> RestartOut
    P --> RestartOut
    T --> RestartOut
    USS --> RestartOut
    VSS --> RestartOut
    PSS --> RestartOut
    TSS --> RestartOut
    
    GX --> SolnFile
    GY --> SolnFile
    GX --> SnapFiles
    GY --> SnapFiles
    
    style Input fill:#E8F5E9
    style Parse fill:#FFF3E0
    style Grid fill:#E3F2FD
    style Init fill:#F3E5F5
    style SolnVars fill:#FFEBEE
    style GridVars fill:#E0F2F1
    style BCVars fill:#FFF8E1
    style Momentum fill:#E1F5FE
    style Pressure fill:#FCE4EC
    style Thermal fill:#F1F8E9
    style SmallScale fill:#EDE7F6
    style Utility fill:#E8EAF6
    style Output fill:#FFE0B2
```

## Data Flow Description

### Input Phase

1. **Parameter Files** → `ReadParam` → Solution variables, flags, tolerances
2. **Boundary Condition File** → `ReadBC` → Boundary condition arrays
3. **Grid File** → `ReadGrid` → Grid coordinates
4. **Restart File** (optional) → `InitialCond` → Initial solution fields

### Processing Phase

1. **Grid Coordinates** → `Metric` → Metric tensors and Jacobians
2. **Solution Variables + Metrics + BCs** → `nAuxMomentum` → Updated velocities
3. **Velocities** → `Filter` → Filtered velocities (optional)
4. **Velocities** → `Divergence` → Divergence field
5. **Divergence + Metrics** → `Ppe` → Pressure field
6. **Pressure + Velocities** → `Project` → Divergence-free velocities
7. **Velocities + Temperature** → `ThermEnergy` → Updated temperature
8. **Pressure + Temperature** → `EqState` → Density field
9. **Velocities + Temperature** → `SmallScale` → Small-scale fields
10. **Small-Scale Fields** → Added to large-scale fields

### Output Phase

1. **Solution Variables** → Various output routines → Output files
2. **Grid Coordinates** → Included in solution output files
3. **Time-Averaged Quantities** → Post-processed and saved separately

## Key Data Transformations

- **Grid → Metrics**: Physical coordinates → Computational metrics
- **Velocities → Auxiliary Velocities**: Through momentum equations
- **Auxiliary Velocities → Pressure**: Through Poisson equation
- **Auxiliary Velocities + Pressure → Divergence-Free Velocities**: Through projection
- **Large-Scale + Small-Scale → Complete Solution**: Through ATD decomposition

