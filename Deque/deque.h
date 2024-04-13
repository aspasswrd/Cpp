#include <iostream>
#include <vector>
#include <stdexcept>

template<typename Type>
class Deque {
 public:
  template<typename IterType>
  class base_iterator;

  friend base_iterator<Type>;

  using iterator = base_iterator<Type>;
  using const_iterator = base_iterator<const Type>;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  iterator begin();
  iterator end();
  const_iterator begin() const { return cbegin(); }
  const_iterator end() const { return cend(); }
  const_iterator cbegin() const;
  const_iterator cend() const;
  reverse_iterator rbegin();
  reverse_iterator rend();
  const_reverse_iterator rbegin() const { return crbegin(); }
  const_reverse_iterator rend() const { return crend(); }
  const_reverse_iterator crbegin() const;
  const_reverse_iterator crend() const;

  Deque();
  explicit Deque(int size);
  Deque(const Deque<Type>& other);
  Deque(int size, const Type& element);

  Deque<Type>& operator=(const Deque<Type>& other) {
    blocks_ = allocate(other.capacity_);
    capacity_ = other.capacity_;
    begin_ = other.begin_;
    size_ = other.size_;
    copy(other.blocks_, blocks_, begin_, begin_ + size_, begin_);
    return *this;
  }

  size_t size() const { return size_; }

  void push_back(const Type& new_element);
  void push_front(const Type& new_element);
  void pop_back();
  void pop_front();

  template<typename IterType>
  void insert(const base_iterator<IterType>& place, const Type& new_element);

  template<typename IterType>
  void erase(const base_iterator<IterType>& place);

  const Type& operator[](size_t index) const {
    return blocks_[(begin_ + index) / kBlock][(begin_ + index) % kBlock];
  }

  Type& operator[](size_t index) {
    return blocks_[(begin_ + index) / kBlock][(begin_ + index) % kBlock];
  }

  Type& at(size_t pos) {
    try {
      if (pos >= size_) {
        throw std::out_of_range("Выход за пределы контейнера");
      }
    }
    catch (...) {
      throw;
    }
    return (*this)[pos];
  }

  const Type& at(size_t pos) const {
    try {
      if (pos >= size_) {
        throw std::out_of_range("Выход за пределы контейнера");
      }
    }
    catch (...) {
      throw;
    }
    return (*this)[pos];
  }

  size_t capacity() { return capacity_; }

  ~Deque();

 private:
  static const ptrdiff_t kBlock = 256;
  std::vector<Type*> blocks_;
  size_t size_;
  size_t begin_;
  size_t capacity_;

  static void copy(const std::vector<Type*>& blocks_old, std::vector<Type*>& blocks_new, size_t begin, size_t end, size_t place);
  static std::vector<Type*> allocate(size_t capacity_new);
  static void deallocate(std::vector<Type*>& blocks, size_t begin, size_t end);

  void increase_capacity();
  void decrease_capacity();
};

template <typename Type>
std::vector<Type*> Deque<Type>::allocate(size_t capacity_new) {
  std::vector<Type*> blocks_new;
  try {
    for (size_t index = 0; index < (capacity_new + kBlock - 1) / kBlock; ++index)
      blocks_new.push_back(reinterpret_cast<Type*>(new char[kBlock * sizeof(Type)]));
  }
  catch (...) {
    for (auto& block : blocks_new)
      delete[] reinterpret_cast<char*>(block);
    throw;
  }
  return blocks_new;
}

template<typename Type>
void Deque<Type>::increase_capacity() {
  if (capacity_ == 0) {
    blocks_ = allocate(3 * kBlock);
    capacity_ = 3 * kBlock;
    begin_ = kBlock;
    return;
  }
  auto new_blocks = allocate(3 * capacity_);
  capacity_ = new_blocks.size() * kBlock;
  size_t place = capacity_ / 3 / kBlock;
  size_t amount = size_ / kBlock + 1;
  for (size_t index = begin_ / kBlock; index < begin_ / kBlock + amount; ++index) {
    delete[] reinterpret_cast<char*>(new_blocks[place]);
    new_blocks[place++] = blocks_[index];
  }
  blocks_ = new_blocks;
  begin_ = capacity_ / 3 / kBlock * kBlock + begin_ % kBlock;
}

template<typename Type>
void Deque<Type>::decrease_capacity() {
  if (blocks_.size() <= 3) {
    return;
  }
  auto new_blocks = allocate(capacity_ / 3);
  capacity_ = new_blocks.size() * kBlock;
  size_t place = capacity_ / 3 / kBlock;
  size_t amount = size_ / kBlock + 2;
  for (size_t index = begin_ / kBlock; index < begin_ / kBlock + amount; ++index) {
    delete[] reinterpret_cast<char*>(new_blocks[place]);
    new_blocks[place++] = blocks_[index];
  }
  blocks_ = new_blocks;
  begin_ = capacity_ / 3 / kBlock * kBlock + begin_ % kBlock;
}

template <typename Type>
void Deque<Type>::deallocate(std::vector<Type*>& blocks, size_t begin, size_t end) {
  Type* block;
  for (size_t index = begin; index < end; ++index) {
    if (index % kBlock == 0 || index == begin)
      block = blocks[index / kBlock];
    //(block + index % kBlock)->~Type();
  }
  for (auto& item : blocks)
    delete[] reinterpret_cast<char*>(item);
  blocks.resize(0);
}

template <typename Type>
void Deque<Type>::copy(const std::vector<Type*>& blocks_old, std::vector<Type*>& blocks_new, size_t begin, size_t end, size_t place) {
  Type* block_old;
  Type* block_new;
  for(size_t index = begin; index < end; ++index, ++place) {
    if (index % kBlock == 0 || index == begin)
      block_old = blocks_old[index / kBlock];

    if (place % kBlock == 0 || index == begin)
      block_new = blocks_new[place / kBlock];

    new (block_new + place % kBlock) Type(block_old[index % kBlock]);
  }
}

template<typename Type>
Deque<Type>::Deque() : blocks_(allocate(0)), size_(0), begin_(0), capacity_(0) {}

template<typename Type>
Deque<Type>::Deque(int size) {
  blocks_ = allocate(3 * size);
  capacity_ = blocks_.size() * kBlock;
  begin_ = capacity_ / 2 - static_cast<size_t>(size) / 2;
  size_ = static_cast<size_t>(size);
  for (int i = 0; i != size; ++i) {
    (*this)[i] = Type();
  }
}

template<typename Type>
Deque<Type>::Deque(const Deque<Type>& other) {
  blocks_ = allocate(other.size_ * 2);
  size_ = other.size_;
  capacity_ = blocks_.size() * kBlock;
  begin_ = (capacity_ - other.size_) / 2;
  copy(other.blocks_, blocks_, other.begin_, other.begin_ + size_, begin_);
}

template<typename Type>
Deque<Type>::Deque(int size, const Type& element) {
  blocks_ = allocate(3 * size);
  capacity_ = blocks_.size() * kBlock;
  size_ = size;
  begin_ = capacity_ / 2 - static_cast<size_t>(size) / 2;;
  for (int i = 0; i != size; ++i) {
    (*this)[i] = element;
  }
}

template <typename Type>
Deque<Type>::~Deque()
{
  deallocate(blocks_, begin_, begin_ + size_);
}

template<typename Type>
void Deque<Type>::push_back(const Type& new_element) {
  if (begin_ + size_ + 2 == capacity_ || capacity_ == 0) {
    increase_capacity();
  }
  blocks_[(begin_ + size_) / kBlock][(begin_ + size_) % kBlock] = new_element;
  ++size_;
}

template<typename Type>
void Deque<Type>::push_front(const Type& new_element) {
  if (begin_ == 0) {
    increase_capacity();
  }
  --begin_;
  blocks_[begin_ / kBlock][begin_ % kBlock] = new_element;
  ++size_;
}

template<typename Type>
void Deque<Type>::pop_back() {
  if (size_ == 0) {
    return;
  }
  --size_;
  if (27 * size_ < capacity_) {
    decrease_capacity();
  }
}

template<typename Type>
void Deque<Type>::pop_front() {
  if (size_ == 0) {
    return;
  }
  ++begin_;
  --size_;
  if (27 * size_ < capacity_) {
    decrease_capacity();
  }
}

template<typename Type>
Deque<Type>::iterator Deque<Type>::begin() {
  if (size_ == 0) {
    return iterator(nullptr);
  }
  return iterator(&blocks_[begin_ / kBlock], &blocks_[begin_ / kBlock][begin_ % kBlock]);
}

template<typename Type>
Deque<Type>::iterator Deque<Type>::end() {
  if (size_ == 0) {
    return iterator(nullptr);
  }
  return iterator(&blocks_[(begin_ + size_) / kBlock], &blocks_[(begin_ + size_) / kBlock][(begin_ + size_) % kBlock]);
}

template<typename Type>
Deque<Type>::const_iterator Deque<Type>::cbegin() const {
  if (size_ == 0) {
    return const_iterator(nullptr);
  }
  return const_iterator(&blocks_[begin_ / kBlock], &blocks_[begin_ / kBlock][begin_ % kBlock]);
}

template<typename Type>
Deque<Type>::const_iterator Deque<Type>::cend() const {
  if (size_ == 0) {
    return const_iterator(nullptr);
  }
  return const_iterator(&blocks_[(begin_ + size_) / kBlock], &blocks_[(begin_ + size_) / kBlock][(begin_ + size_) % kBlock]);
}

template<typename Type>
Deque<Type>::reverse_iterator Deque<Type>::rend() {
  if (size_ == 0) {
    return reverse_iterator(iterator(nullptr));
  }
  return reverse_iterator(iterator(&blocks_[begin_ / kBlock], &blocks_[begin_ / kBlock][begin_ % kBlock]));
}

template<typename Type>
Deque<Type>::reverse_iterator Deque<Type>::rbegin() {
  if (size_ == 0) {
    return reverse_iterator(iterator(nullptr));
  }
  return reverse_iterator(iterator(&blocks_[(begin_ + size_) / kBlock], &blocks_[(begin_ + size_) / kBlock][(begin_ + size_) % kBlock]));
}

template<typename Type>
Deque<Type>::const_reverse_iterator Deque<Type>::crbegin() const {
  if (size_ == 0) {
    return const_reverse_iterator(const_iterator(nullptr));
  }
  return const_reverse_iterator(const_iterator(&blocks_[(begin_ + size_) / kBlock], &blocks_[(begin_ + size_) / kBlock][(begin_ + size_) % kBlock]));
}

template<typename Type>
Deque<Type>::const_reverse_iterator Deque<Type>::crend() const {
  if (size_ == 0) {
    return const_reverse_iterator(const_iterator(nullptr));
  }
  return const_reverse_iterator(const_iterator(&blocks_[begin_ / kBlock], &blocks_[begin_ / kBlock][begin_ % kBlock]));
}

template <typename Type>
template<typename IterType>
void Deque<Type>::insert(const base_iterator<IterType>& place, const Type& new_element) {
  if (place == end())
  {
    push_back(new_element);
    return;
  }
  if (place == begin())
  {
    push_front(new_element);
    return;
  }
  std::vector<Type*> blocks_new;
  try
  {
    blocks_new = allocate(capacity_ + 1);
    size_t index = place - begin();
    copy(blocks_, blocks_new, begin_, begin_ + index, begin_);
    new (blocks_new[(begin_ + index) / kBlock] + (begin_ + index) % kBlock) Type(new_element);
    copy(blocks_, blocks_new, begin_ + index, begin_ + size_, begin_ + index + 1);
  }
  catch(...)
  {
    for (auto& block : blocks_new)
      delete[] reinterpret_cast<char*>(block);
    return;
  }
  deallocate(blocks_, begin_, begin_ + size_);
  blocks_ = blocks_new;
  ++size_;
  ++capacity_;
}

template<typename Type>
template<typename IterType>
void Deque<Type>::erase(const Deque::base_iterator<IterType>& place) {
  if (place == begin()) {
    pop_front();
    return;
  }
  if (place == end()) {
    pop_back();
    return;
  }
  std::vector<Type*> blocks_new;
  try {
    blocks_new = allocate(capacity_ - 1);
    size_t index = place - begin();
    copy(blocks_, blocks_new, begin_, begin_ + index - 1, begin_);
    copy(blocks_, blocks_new, begin_ + index, begin_ + size_, begin_ + index - 1);
  }
  catch(...) {
    for (auto& block : blocks_new)
      delete[] reinterpret_cast<char*>(block);
    return;
  }
  deallocate(blocks_, begin_, begin_ + size_);
  blocks_ = blocks_new;
  --size_;
  --capacity_;
}


template<typename Type>
template<typename IterType>
class Deque<Type>::base_iterator {
 public:
  friend Deque<Type>;
  using value_type = IterType;
  using difference_type = ptrdiff_t;
  using pointer = IterType*;
  using reference = IterType&;
  using iterator_category = std::random_access_iterator_tag;

  base_iterator() = default;
  base_iterator(IterType* const* block, IterType* elem) : block_(block), elem_(elem), elem_first_(*block) {};
  //base_iterator(const base_iterator<IterType>& other) : block_(other.block_), elem_(other.elem_), elem_first_(other.elem_first_) {}
  explicit base_iterator(nullptr_t) : block_(nullptr), elem_(nullptr), elem_first_(nullptr) {}

  base_iterator<IterType>& operator+=(difference_type diff) {
    if (elem_ == nullptr) {
      return *this;
    }
    if (diff < 0) {
      return (*this -= -diff);
    }
    if (diff == 0) {
      return *this;
    }
    difference_type new_index = diff + (elem_ - elem_first_);
    if (new_index - kBlock >= 0) {
      block_ += new_index / kBlock;
      elem_first_ = *block_;
    }
    elem_ = elem_first_ + new_index % kBlock;
    return *this;
  }

  base_iterator<IterType> operator+(difference_type diff) const {
    return base_iterator<IterType>(*this) += diff;
  }

  base_iterator<IterType> operator-(difference_type diff) const {
    return base_iterator<IterType>(*this) -= diff;
  }

  difference_type operator-(const base_iterator<IterType>& other) const {
    if (block_ == other.block_) {
      return elem_ - other.elem_;
    }
    long difference = block_ - other.block_;
    difference -= 1;
    long to_back = kBlock - (other.elem_ - other.elem_first_);
    long to_start = elem_ - elem_first_;
    return (difference * kBlock + to_back + to_start);
  }

  base_iterator<IterType>& operator-=(difference_type diff) {
    if (elem_ == nullptr) {
      return *this;
    }
    if (diff < 0) {
      return (*this += -diff);
    }
    if (diff == 0) {
      return *this;
    }
    difference_type index_in_block = (elem_ - elem_first_);
    while (diff > index_in_block) {
      block_--;
      elem_first_ = *block_;
      index_in_block += kBlock;
    }
    elem_ = elem_first_ + (index_in_block - diff);
    return *this;
  }


  base_iterator<IterType>& operator++() { return *this += 1; }
  base_iterator<IterType> operator++(int) { return *this += 1; }
  IterType* operator->() { return elem_; }
  const IterType* operator->() const { return elem_; }
  base_iterator<IterType>& operator--() { return (*this -= 1); }
  const IterType& operator*() const { return *elem_; }
  IterType& operator*() { return *elem_; }

  bool operator<(const base_iterator<IterType>& other) const {
    if (block_ == other.block_) {
      return (elem_ - elem_first_) < (other.elem_ - other.elem_first_);
    }
    return block_ < other.block_;
  }
  bool operator>(const base_iterator<IterType>& other) const { return other < *this; }
  bool operator==(const base_iterator<IterType>& other) const { return elem_ == other.elem_; }
  bool operator>=(const base_iterator<IterType>& other) const { return elem_ == other.elem_ || *this > other; }
  bool operator<=(const base_iterator<IterType>& other) const { return elem_ == other.elem_ || *this < other; }
  bool operator!=(const base_iterator<IterType>& other) const { return elem_ != other.elem_; }
  //base_iterator<IterType>& operator=(const Type* ptr) { elem_ = ptr; }

  operator base_iterator<const IterType>() const {
    return base_iterator<const IterType>(*this);
  }

 private:
  IterType* const* block_;
  IterType* elem_;
  IterType* elem_first_;
  static const ptrdiff_t kBlock = Deque<Type>::kBlock;
};

