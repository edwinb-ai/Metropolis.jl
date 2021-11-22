abstract type GeneralCache end

mutable struct CellCache{C,A} <: GeneralCache
    cell::C
    aux::A
end

function _build_cache(cell::CellListMap.CellList)
    aux = CellListMap.AuxThreaded(cell)

    return CellCache(cell, aux)
end
