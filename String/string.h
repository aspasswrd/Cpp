#include <iostream>
#include <cstring>

class String {
 public:
  String();
  ~String() { delete[] arr_; };
  String(const char*);
  String(char* arr, size_t size, size_t capacity);
  String(int length, char c);
  String(const String&);

  void shrink_to_fit();
  void clear();

  String& operator= (const String& str);
  String& operator= (const char* s);
  bool operator== (const String& str) const;
  bool operator!= (const String& str) const;
  bool operator< (const String& str) const;
  bool operator> (const String& str) const;
  bool operator<=(const String& str) const {
    return (*this) < str || (*this) == str;
  }
  bool operator>= (const String& str) const {
    return (*this) > str || (*this) == str;
  }
  char& operator[] (int index);
  char& operator[] (size_t index);
  const char& operator[] (int index) const;
  const char& operator[] (size_t index) const;

  const size_t& length() const { return size_; }
  const size_t& size() const { return size_; }
  const size_t& capacity() const { return capacity_; }
  bool empty() const { return size_ == 0; }

  void push_back(char);
  void pop_back();

  char* begin() const { return arr_; }
  char* end() const { return &arr_[size_ + 1]; }
  char* data() const { return arr_; }

  char& front() { return arr_[0]; }
  const char& front() const { return arr_[0]; }
  char& back() { return arr_[size_ - 1]; }
  const char& back() const { return arr_[size_ - 1]; }

  String operator+ (const String& str) const;
  String operator+ (char) const;

  String& operator+= (const String& str);
  String& operator+= (char c);

  size_t find(const String& substring) const;
  size_t rfind(const String& substring) const;
  String substr(size_t pos, size_t len) const;

 private:
  static size_t nearest_deg_of_2(size_t num) {
    while (!(num > 0 && (num & (num - 1)) == 0)) {
      num++;
    }
    return num;
  }

  static void IncreaseCapacity(char*& arr, size_t& capacity) {
    char* resize = new char [(capacity + 1) * 2];
    std::copy(arr, &arr[capacity], resize);
    delete[] arr;
    arr = resize;
    capacity = (capacity + 1) * 2 - 1;
  }

  char* arr_;
  size_t size_;
  size_t capacity_;
};

String::String() {
  arr_ = new char [8];
  arr_[0] = '\0';
  size_ = 0;
  capacity_ = 7;
}

String::String(const char* str) {
  size_t len = strlen(str);
  size_t deg = nearest_deg_of_2(len + 1);
  arr_ = new char [deg];
  std::copy(str, &str[len] ,arr_);
  arr_[len] = '\0';
  size_ = len;
  capacity_ = deg - 1;
}

String::String(char* arr, size_t size, size_t capacity) : arr_(arr)
                                                        , size_(size)
                                                        , capacity_(capacity) {}

String::String(int length, char c) {
  size_t deg = nearest_deg_of_2(length);
  arr_ = new char [deg];
  arr_[length] = '\0';
  std::fill(arr_, &arr_[length], c);
  size_ = length;
  capacity_ = deg - 1;
}

String::String(const String& str) {
  arr_ = new char [str.capacity_ + 1];
  std::copy(str.begin(), str.end(), arr_);
  size_ = str.size_;
  capacity_ = str.capacity_;
}

void String::shrink_to_fit() {
  char* resize = new char [size_ + 1];
  std::copy(this->begin(), this->end(), resize);
  delete[] arr_;
  arr_ = resize;
  capacity_ = size_;
}

void String::clear() {
  for (size_t i = 0; i != size_; ++i) {
    arr_[i] = '\0';
  }
  size_ = 0;
}

String& String::operator=(const String& str) {
  if (str == *this || this == &str) { return *this; }
  delete[] arr_;
  arr_ = new char [str.capacity_ + 1];
  std::copy(str.begin(), str.end(), arr_);
  capacity_ = str.capacity_;
  size_ = str.size_;
  return *this;
}

String& String::operator=(const char* s) {
  delete[] arr_;
  size_t len = nearest_deg_of_2(strlen(s) + 1);
  arr_ = new char [len];
  std::copy(s, &s[strlen(s)], arr_);
  size_ = strlen(s);
  capacity_ = len - 1;
  return *this;
}

bool String::operator==(const String &str) const {
  if (str.size_ != size_) { return false; }
  for (size_t i = 0; i != str.size_; ++i) {
    if (arr_[i] != str.arr_[i]) { return false; }
  }
  return true;
}

bool operator==(const char* str1, const String& str2) {
  if (strlen(str1) != str2.size()) { return false; }
  for (size_t i = 0; i != strlen(str1); ++i) {
    if (str1[i] != str2[i]) { return false; }
  }
  return true;
}

bool String::operator!=(const String& str) const {
  return !(str == (*this));
}

bool operator!=(const char* str1, const String& str2) {
  return !(str1 == str2);
}

bool String::operator<(const String &str) const {
  return strcmp(this->data(), str.data()) < 0;
}

bool operator<(const char* str1, const String& str2) {
  return strcmp(str1, str2.data()) < 0;
}

bool String::operator>(const String &str) const {
  return strcmp(this->data(), str.data()) > 0;
}

bool operator>(const char* str1, const String& str2) {
  return strcmp(str1, str2.data()) > 0;
}

bool operator<=(const char* str1, const String& str2) {
  return str1 == str2 || str1 < str2;
}

bool operator>=(const char* str1, const String& str2) {
  return str1 == str2 || str1 > str2;
}

char& String::operator[](int index) { return arr_[index]; }
char& String::operator[](size_t index) { return arr_[index]; }

const char& String::operator[](int index) const {
  return arr_[index];
}

const char& String::operator[](size_t index) const {
  return arr_[index];
}

void String::push_back(const char c) {
  if (size_ == capacity_) {
    IncreaseCapacity(arr_, capacity_);
  }
  arr_[size_] = c;
  size_++;
  arr_[size_] = '\0';
}

void String::pop_back() {
  if (size_ == 0) { return; }
  size_--;
  arr_[size_] = '\0';
}

String String::operator+(const String &str) const {
  String answer;
  delete[] answer.arr_;
  auto len1 = size_;
  auto len2 = strlen(str.arr_);
  answer.capacity_ = nearest_deg_of_2(len1 + len2 + 1) - 1;
  answer.arr_ = new char [answer.capacity_ + 1];
  answer.size_ = len1 + len2;
  std::copy(this->begin(), this->end(), answer.arr_);
  std::copy(str.begin(), str.end(), &answer.arr_[len1]);
  return answer;
}

String String::operator+(char c) const {
  String answer(*this);
  answer.push_back(c);
  return answer;
}

String operator+(char c, const String& str) {
  String temp(1, c);
  return temp + str;
}

String operator+(const char* c, const String& str) {
  String answer(c);
  answer += str;
  return answer;
}

String& String::operator+=(char c) {
  this->push_back(c);
  return *this;
}

String& String::operator+=(const String &str) {
  size_t len1 = size_;
  size_t len2 = str.size_;
  size_t new_capacity = nearest_deg_of_2(len1 + len2 + 1);
  char* resize = new char [new_capacity];
  std::copy(arr_, &arr_[len1], resize);
  delete[] arr_;
  arr_ = resize;
  capacity_ = new_capacity - 1;
  std::copy(str.begin(), str.end(), &arr_[len1]);
  size_ = len1 + len2;
  return *this;
}

String String::substr(size_t pos, size_t len) const {
  if (pos + len > size_) { len = size_ - pos; }
  char* answer = new char[len + 1];
  std::copy(&arr_[pos], &arr_[pos + len], answer);
  arr_[len] = '\0';
  return {answer, len, len};
}

size_t String::find(const String& substring) const {
  bool flag = false;
  for (size_t i = 0; i != size_ - substring.size() + 1; ++i) {
    for (size_t j = 0; j != substring.size(); ++j) {
      flag = true;
      if (arr_[i + j] != substring[j]) {
        flag = false;
        break;
      }
    }
    if (flag) {
      return i;
    }
  }
  return size_;
}

size_t String::rfind(const String& substring) const {
  bool flag = false;
  size_t answer = size_;
  for (size_t i = 0; i != size_ - substring.size() + 1; ++i) {
    for (size_t j = 0; j != substring.size(); ++j) {
      flag = true;
      if (arr_[i + j] != substring[j]) {
        flag = false;
        break;
      }
    }
    if (flag) {
      answer = i;
    }
  }
  return answer;
}

std::istream& operator>>(std::istream& in, String& str) {
  str.clear();
  int c = in.get();
  while (c != '\n' && c != ' ' &&  !in.eof()) {
    str.push_back(static_cast<char>(c));
    c = in.get();
  }
  return in;
}

std::ostream& operator<<(std::ostream& out, const String& str) {
  out << str.data();
  return out;
}


