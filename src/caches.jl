"""
    GeneralCache

Supertype for all kinds of cache types that are needed throughout the simulation.
"""
abstract type GeneralCache end

"""
    CellCache{C,A} <: GeneralCache

A special type of cache to handle `CellList` and auxiliary object.

# Fields
- `cell::C`: stores the cell list object, normally of type `CellListMap.CellList`.
- `aux::A`: stores the auxiliary lists for cell lists, useful when having to continously
    update the `cell` field.
"""
mutable struct CellCache{A} <: GeneralCache
    cell::CellListMap.CellList
    aux::A
end

function _build_cache(cell::CellListMap.CellList)
    aux = CellListMap.AuxThreaded(cell)

    return CellCache(cell, aux)
end
