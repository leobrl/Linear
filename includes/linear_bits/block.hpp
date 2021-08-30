#pragma once

#include<vector>
#include<memory>

#include <stdlib.h>
#include <new>
#include <limits>


namespace linear{
    
    template<typename T> class AlignAllocator;
    template <typename T> class Block;
    template <typename T> std::ostream& operator<<( std::ostream&, const Block<T>& );

    template<typename T>
    class Block{
        public:

            using raw_memory    = std::vector<T, AlignAllocator<T> >;
            using raw_memory_it = typename std::vector<T, AlignAllocator<T> >::iterator;

        protected:

            raw_memory buffer;

        public:

            inline                  Block()  = default;
            inline                  ~Block() = default;

            inline explicit         Block(const size_t);
            inline explicit         Block(const raw_memory& );
            inline explicit         Block(const std::vector<T>& );

                                    Block(const Block&) = default;
                                    Block(Block&&)      = default;

            Block&                  operator= (const Block&) = default;
            Block&                  operator= (Block&&)      = default;    

            void                    fill(const T&);
            void                    fill(const raw_memory&);
            void                    fill(const std::vector<T>&);
            
            inline const T&         operator[] (size_t idx) const   { return buffer[idx]; }
            inline T&               operator[] (size_t idx)         { return buffer[idx]; }

            inline size_t           size () const   { return buffer.size(); }

            inline raw_memory_it    begin()     {return buffer.begin();}
            inline raw_memory_it    end()       {return buffer.end();}

            inline T*               pfront() {return buffer.data();}

            // Friends
            friend std::ostream& operator<< <T> (std::ostream&, const Block<T>&);
    };

    template <typename T>
    class AlignAllocator
    {
        public:
            typedef T value_type;
            
            AlignAllocator () = default;

            template<typename U> AlignAllocator (const AlignAllocator<U> &) noexcept {}
        
//            template <typename U> struct rebind { using other = AlignAllocator<U, T_num_els>; };

            template<class U> bool operator==(const AlignAllocator<U>&) const noexcept
            {
                return true;
            }
        
            template<class U> bool operator!=(const AlignAllocator<U>&) const noexcept
            {
                return false;
            }

            [[nodiscard]] T* allocate(std::size_t n) {
                size_t alignment = 16;

                if (n > std::numeric_limits<std::size_t>::max() / sizeof(T))
                throw std::bad_alloc();
                
                if (auto p = static_cast<T*>(std::aligned_alloc(alignment, n*sizeof(T))))
                return p;
        
                throw std::bad_alloc();
            }
            
            void deallocate(T* const p, std::size_t) const noexcept { std::free(p); }
    };
}
