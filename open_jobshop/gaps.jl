type Gap
    start::Int64
    size::Int64
end

typealias GapList Vector{Gap}

# Index of first slot *after* the gap
fin(gap::Gap) = gap.start + gap.size

# Useful for functions like max, sort:
isless(gap1::Gap, gap2::Gap) = isless(gap1.start, gap2.start)


#
# NOT USED
#
# function intersect(gaps1::GapList, gaps2::GapList)

#     gaps = Gap[]
    
#     for gap1 in gaps1
#         for gap2 in gaps2
#             # Case 1
#             # 1 -------oooo--
#             # 2 -ooo---------
#             if gap1.start >= fin(gap2)
#                 continue
#             # Case 2
#             # 1 -oooo--------
#             # 2 -------ooo---
#             elseif gap2.start >= fin(gap1)
#                 break # gaps are sorted, so no need to continue
#             else

#                 start = max([gap1.start, gap2.start])
#                 after = min([fin(gap1), fin(gap2)])
#                 size = after - start

#                 push(gaps, Gap(start, size))
#             end

#         end # loop over gaps2
#     end # loop over gaps1

#     sort!(gaps)
#     return gaps
# end


# Subtract a gap from a gap list, meaning gaps will be closed, e.g.
#----ooooo----
# -
#-----oo------
# =
#----o--oo----
function subtract(gaps::GapList, gap::Gap)

    new_gaps = Gap[]
    
    for i = 1:length(gaps) 
        if gaps[i].start >= fin(gap) || gap.start >= fin(gaps[i])
            push(new_gaps, gaps[i])
        else
            #del(gaps, i)
            if gap.start > gaps[i].start
                left_gap_size = gap.start - gaps[i].start 
                push(new_gaps, Gap(gaps[i].start, left_gap_size))
            end

            if fin(gap) < fin(gaps[i])
                right_gap_size = fin(gaps[i]) - fin(gap)
                push(new_gaps, Gap(fin(gap), right_gap_size)) 
            end
        end

    end # loop over gaps

    return new_gaps
end



# close the first available gap
function find_gap(list::GapList, size)

    for i = 1:length(list)
        gap = list[i]
        if gap.size >= size
            return gap.start
            break
        end
    end

    return 0
end


# function test()

#     for i  = 9:12
#         l = [Gap(1,3), Gap(5,2), Gap(9,4)]
#         l2 = subtract(l, Gap(i,1))
#         println(l2)
#     end

# end

# test()

