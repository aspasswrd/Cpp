#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <random>

namespace {
    const double precision = 1e-6;
    const double PI = 3.14159265358979324;

    bool equal_double(double d1, double d2) {
        return std::abs(d1 - d2) < precision;
    }

    bool less(double d1, double d2) {
        return d1 - d2 < -precision;
    }

    bool lessOrEqual(double d1, double d2) {
        return less(d1, d2) || equal_double(d1, d2);
    }
}

struct Point {
  double x = 0;
  double y = 0;

  Point() = default;

  Point(double qx, double qy) {
    x = qx;
    y = qy;
  }

  bool operator==(const Point& point) const {
    if (equal_double(point.x, x) && equal_double(point.y, y)) {
      return true;
    }
    return false;
  }

  bool operator!=(const Point& point) const {
    return !(*this == point);
  }

  void rotate(const Point& center, double phi) {
    phi = (PI / 180) * phi;

    x -= center.x;
    y -= center.y;

    double X = x * cos(phi) + y * sin(phi);
    double Y = (-1) * x * sin(phi) + y * cos(phi);

    x = X + center.x;
    y = Y + center.y;
  }
};

struct Vec {
  double dx = 0;
  double dy = 0;

  Vec() = default;

  Vec(const Point& p1, const Point& p2) {
    dx = p2.x - p1.x;
    dy = p2.y - p1.y;
  }

  Vec(double x, double y) : dx(x), dy(y) {}

  double Length() const {
    return std::sqrt(dx * dx + dy * dy);
  }

  double LengthSqr() const {
    return dx * dx + dy * dy;
  }

  Vec& operator/= (double value) {
    dx /= value;
    dy /= value;
    return *this;
  }

  Vec& operator*= (double value) {
    dx *= value;
    dy *= value;
    return *this;
  }

  void Normalize() {
    double coef = Length();
    dx /= coef;
    dy /= coef;
  }

  Vec ort() const {
    return {this->dy * (-1), this->dx};
  }
};

Vec operator*(const Vec& vec, double val) {
  Vec answer = Vec(vec.dx * val, vec.dy * val);
  return answer;
}

Point operator+ (const Point& p, const Vec& vec) {
  return Point{p.x + vec.dx, p.y + vec.dy};
}

Point operator- (const Point& p, const Vec& vec) {
  return Point{p.x - vec.dx, p.y - vec.dy};
}

Vec operator+ (const Vec& v1, const Vec& v2) {
  return {v1.dx + v2.dx, v1.dy + v2.dy};
}

bool SignVecProduct(const Vec& vec1, const Vec& vec2) {
  double product = vec1.dx * vec2.dy - vec1.dy * vec2.dx;
  return product > 0;
}

double ScalarProduct(const Vec& vec1, const Vec& vec2) {
  return vec1.dx * vec2.dx + vec1.dy * vec2.dy;
}

bool isCollinear(const Vec& vec1, const Vec& vec2) {
  return std::abs(vec1.dx * vec2.dy - vec2.dx * vec1.dy) < precision;
}

class Line {
 public:
  Line() = default;

  Line(const Point& p1, const Point& p2) {
    Vec GuideVec = Vec(p1, p2);
    b_ = -(GuideVec.dx);
    if (equal_double(b_, 0.0)) {
      b_ = +0.0;
    }
    a_ = GuideVec.dy;
    if (equal_double(a_, 0.0)) {
      a_ = +0.0;
    }
    c_ = -(a_ * p1.x + b_ * p1.y);
    if (equal_double(c_, 0.0)) {
      c_ = +0.0;
    }
  }

  Line(double k, double b) {
    //y = kx + b
    a_ = k;
    if (equal_double(a_, 0.0)) {
      a_ = +0.0;
    }
    b_ = -1;
    if (equal_double(b_, 0.0)) {
      b_ = +0.0;
    }
    c_ = b;
    if (equal_double(c_, 0.0)) {
      c_ = +0.0;
    }
  }

  Line(const Point& p, double k) {
    a_ = k;
    b_ = -1;
    c_ = p.y - k * p.x;
  }

  Line(const Point& p, const Vec& vec){
    b_ = -vec.dx;
    if (equal_double(b_, 0.0)) {
      b_ = +0.0;
    }
    a_ = vec.dy;
    if (equal_double(a_, 0.0)) {
      a_ = +0.0;
    }
    c_ = -(a_ * p.x + b_ * p.y);
    if (equal_double(c_, 0.0)) {
      c_ = +0.0;
    }
  }

  Vec GuideVec() const { return Vec{b_ * (-1), a_}; }

  bool ContainsPoint(const Point& p) const {
    return equal_double(a_ * p.x + b_ * p.y + c_, 0.0);
  }

  bool operator==(const Line& line) const {
    Vec vec1 = this->GuideVec();
    Vec vec2 = line.GuideVec();
    if (!isCollinear(vec1, vec2)) {
      return false;
    }
    return a_ * line.c_ - line.a_ * c_ < precision;
  }

  double GetA() const { return a_; }
  double GetB() const { return b_; }
  double GetC() const { return c_; }

 private:
  // Ax + By + C = 0
  double a_;
  double b_;
  double c_;
};

double DistanceToLine(const Point& p, const Line& l) {
  double A = l.GetA();
  double B = l.GetB();
  double C = l.GetC();

  double answer = (A * p.x + B * p.y + C) / std::sqrt(A * A + B * B);
  return std::abs(answer);
}

Point Intersection(const Line& l1, const Line& l2) {
  double A1 = l1.GetA();
  double B1 = l1.GetB();
  double C1 = l1.GetC();

  double A2 = l2.GetA();
  double B2 = l2.GetB();
  double C2 = l2.GetC();

  double Det = A1 * B2 - A2 * B1;

  if (equal_double(Det, 0)) {
    throw std::logic_error("Линии параллельны");
  } else {
    double DetX = B1 * C2 - B2 * C1;
    double DetY = C1 * A2 - C2 * A1;
    Point answer(DetX / Det, DetY / Det);
    if (equal_double(answer.x, 0.0)) {
      answer.x = +0.0;
    }
    if (equal_double(answer.y, 0.0)) {
      answer.y = +0.0;
    }
    return answer;
  }
}

void p_reflect(Point& p, const Point& center) {
  Vec vec(p, center);
  vec *= 2;
  p = p + vec;
}

void p_reflect(Point& p, const Line& axis) {
  Vec ort = axis.GuideVec().ort();
  ort.Normalize();
  double coef = DistanceToLine(p, axis);
  ort *= coef;
  if (!axis.ContainsPoint(p + ort)) {
    ort *= (-1);
  }
  ort *= 2;
  p = p + ort;
}

void p_scale(Point& p, const Point& center, double coefficient) {
  Vec vec(center, p);
  vec *= coefficient;
  p = center + vec;
}

class Shape {
 public:
  virtual ~Shape() = default;

  virtual double perimeter() const = 0;
  virtual double area() const = 0;

  virtual bool operator==(const Shape& another) const = 0;
  virtual bool operator!=(const Shape& anothe) const = 0;
  virtual bool isCongruentTo(const Shape& another) const = 0;
  virtual bool isSimilarTo(const Shape& another) const = 0;
  virtual bool containsPoint(const Point& point) const = 0;

  virtual void rotate(const Point& center, double angle) = 0;
  virtual void reflect(const Point& center) = 0;
  virtual void reflect(const Line& axis) = 0;
  virtual void scale(const Point& center, double coefficient) = 0;
};

class Polygon : public Shape {
 public:
  Polygon() = default;

  Polygon(const std::vector<Point>& vert) {
    vertices_ = std::vector<Point>(vert.begin(), vert.end());
  }

  template<class ... Points>
  Polygon(Points&& ... points) : vertices_{std::forward<Points>(points)...} {}

  template<>
  Polygon(Polygon& another) {
    std::copy(another.vertices_.begin(), another.vertices_.end(), std::back_inserter(vertices_));
  }

  size_t verticesCount() const { return vertices_.size(); }

  virtual bool isConvex() {
    if (verticesCount() <= 3) {
      return true;
    }

    Vec vec1(vertices_[verticesCount() - 1], vertices_[0]);
    Vec vec2(vertices_[0], vertices_[1]);
    bool sign = SignVecProduct(vec1, vec2);
    bool new_sign;

    for (size_t i = 0; i != verticesCount() - 1; ++i) {
      vec1 = Vec(vertices_[i], vertices_[i + 1]);
      vec2 = Vec(vertices_[i + 1], vertices_[i + 2]);
      new_sign = SignVecProduct(vec1, vec2);
      if (sign != new_sign) {
        return false;
      }
    }
    return true;
  }

  std::vector<Point> getVertices() { return vertices_; }

  double perimeter() const override {
    double answer = Vec(vertices_[verticesCount() - 1], vertices_[0]).Length();
    for (size_t i = 0; i != verticesCount() - 1; ++i) {
      answer += Vec(vertices_[i], vertices_[i + 1]).Length();
    }
    return answer;
  }

  double area() const override {
    if (vertices_.size() <= 2) {
      return 0;
    }
    double answer = 0;
    for (size_t i = 0; i != verticesCount(); ++i) {
      Point p1 = vertices_[i];
      Point p2 = vertices_[(i + 1) % verticesCount()];
      answer += (p2.x - p1.x) * (p2.y + p1.y) / 2;
    }
    if (answer < 0) {
      answer *= -1;
    }
    return answer;
  }

  bool containsPoint(const Point& p) const override {
    std::random_device random_device;
    std::mt19937 gen(random_device());
    std::uniform_int_distribution<> dist(-100000, 100000);
    Line l(p, dist(gen));
    std::vector<Point> IntersectionPoints;

    Line segment(vertices_[0], vertices_[verticesCount() - 1]);
    double maxX = std::max(vertices_[0].x, vertices_[verticesCount() - 1].x);
    double maxY = std::max(vertices_[0].y, vertices_[verticesCount() - 1].y);

    double minX = std::min(vertices_[0].x, vertices_[verticesCount() - 1].x);
    double minY = std::max(vertices_[0].y, vertices_[verticesCount() - 1].y);

    Point InterI = Intersection(l, segment);
    if (lessOrEqual(InterI.x, maxX) && lessOrEqual(InterI.y, maxY) &&
        lessOrEqual(minX, InterI.x) && lessOrEqual(minY, InterI.y) && lessOrEqual(p.x, InterI.x)) {
      IntersectionPoints.push_back(InterI);
    }

    for (size_t i = 0; i != verticesCount() - 1; ++i) {
      segment = Line(vertices_[i], vertices_[i + 1]);

      maxX = std::max(vertices_[i].x, vertices_[i + 1].x);
      maxY = std::max(vertices_[i].y, vertices_[i + 1].y);

      minX = std::min(vertices_[i].x, vertices_[i + 1].x);
      minY = std::min(vertices_[i].y, vertices_[i + 1].y);

      if (isCollinear(l.GuideVec(),segment.GuideVec())) {
        if (equal_double(p.y, vertices_[i].y) && less(p.x, maxX) && less(minX, p.x)) {
          return true;
        }
        continue;
      }
      InterI = Intersection(l, segment);
      if (lessOrEqual(InterI.x, maxX) && lessOrEqual(InterI.y, maxY) &&
          lessOrEqual(minX, InterI.x) && lessOrEqual(minY, InterI.y) && lessOrEqual(p.x, InterI.x)) {
        IntersectionPoints.push_back(InterI);
      }
    }

    IntersectionPoints.erase(
        std::unique(IntersectionPoints.begin(), IntersectionPoints.end()),
        IntersectionPoints.end());
    if (IntersectionPoints.size() == 1) {
      return true;
    }
    return false;
  }

  bool isCongruentTo(const Shape& another) const override {
    const auto* ptr = dynamic_cast<const Polygon*>(&another);
    if (ptr == nullptr) {
      return false;
    }

    if (!isSimilarTo(another)) {
      return false;
    }

    std::vector<double> Lenghts1 = GetLenghts(vertices_);
    std::vector<double> Lenghts2 = GetLenghts(ptr->vertices_);

    return CheckIsomorphic(Lenghts1, Lenghts2);
  }

  bool isSimilarTo(const Shape& another) const override {
    const auto* ptr = dynamic_cast<const Polygon*>(&another);
    if (ptr == nullptr) {
      return false;
    }
    if (this->verticesCount() != ptr->verticesCount()) {
      return false;
    }
    std::vector<double> Cosines1 = GetCosines(vertices_);
    std::vector<double> Cosines2 = GetCosines(ptr->vertices_);

    return CheckIsomorphic(Cosines1, Cosines2);
  }

  bool operator==(const Shape& another) const override {
    const auto* ptr = dynamic_cast<const Polygon*>(&another);
    if (ptr == nullptr) {
      return false;
    }
    return CheckIsomorphic(vertices_, ptr->vertices_);
  }

  bool operator!=(const Shape& another) const override {
    return !(*this == another);
  }

  void rotate(const Point& center, double angle) override {
    for (auto& vert : vertices_) {
      vert.rotate(center, angle);
    }
  }

  void reflect(const Point& center) override {
    for (auto& vert : vertices_) {
      p_reflect(vert, center);
    }
  }

  void reflect(const Line& axis) override {
    for (auto& vert : vertices_) {
      p_reflect(vert, axis);
    }
  }

  void scale(const Point& center, double coefficient) override {
    for (auto& vert : vertices_) {
      p_scale(vert, center, coefficient);
    }
  }

 protected:
  std::vector<Point> vertices_;

 private:
  static std::vector<double> GetLenghts(const std::vector<Point>& vec) {
    std::vector<double> answer(vec.size());
    Vec veci(vec[vec.size() - 1], vec[0]);
    answer[vec.size() - 1] = veci.LengthSqr();

    for (size_t i = 0; i != vec.size() - 1; ++i) {
      veci = Vec(vec[i], vec[i + 1]);
      answer[i] = veci.LengthSqr();
    }

    return answer;
  }

  static std::vector<double> GetCosines(const std::vector<Point>& vec) {
    if (vec.size() == 2) {
      return {1, 1};
    }
    std::vector<double> answer(vec.size());

    Vec veci1(vec[0], vec[vec.size() - 1]);
    Vec veci2(vec[0], vec[1]);
    answer[0] = ScalarProduct(veci1, veci2) / (veci1.Length() * veci2.Length());

    for (size_t i = 1; i != vec.size() - 1; ++i) {
      veci1 = Vec(vec[i], vec[i - 1]);
      veci2 = Vec(vec[i], vec[i + 1]);

      answer[i] = ScalarProduct(veci1, veci2) / (veci1.Length() * veci2.Length());
    }

    veci1 = Vec(vec[vec.size() - 1], vec[0]);
    veci2 = Vec(vec[vec.size() - 1], vec[vec.size() - 2]);

    answer[vec.size() - 1] = ScalarProduct(veci1, veci2) / (veci1.Length() * veci2.Length());

    return answer;
  }

  static bool equal(double d1, double d2) {
    return equal_double(d1, d2);
  }

  static bool equal(Point p1, Point p2) {
    return p1 == p2;
  }

  template<class T>
  bool CheckIsomorphic(const std::vector<T>& vec1, const std::vector<T>& vec2) const {
    return CheckCyclicShift(vec1, vec2) || CheckCyclicShift(std::vector<T>(vec1.rbegin(), vec1.rend()), vec2);
  }

  template<class T>
  bool CheckCyclicShift(const std::vector<T>& vec1, const std::vector<T>& vec2) const {
    if (vec1.size() != vec2.size()) {
      return false;
    }
    std::vector<T> concatenate(2 * vec1.size());
    std::copy(&vec1[0], &vec1[vec1.size()], concatenate.begin());
    std::copy(&vec1[0], &vec1[vec1.size()], &concatenate[vec1.size()]);

    bool flag = false;
    for (size_t i = 0; i != vec1.size() + 1; ++i) {
      for (size_t j = 0; j != vec2.size(); ++j) {
        flag = true;
        if (!equal(concatenate[i + j], vec2[j])) {
          flag = false;
          break;
        }
      }
      if (flag) {
        return true;
      }
    }
    return false;
  }
};

class Ellipse : public Shape {
 public:
  Ellipse(const Point& p1, const Point& p2, double sum) {
    focus_1_ = p1;
    focus_2_ = p2;
    a_ = sum / 2;
    Vec focuses = Vec(p1, p2);
    c_ = focuses.Length() / 2;
    b_ = sqrt(a_ * a_ - c_ * c_);
  }

  std::pair<Point, Point> focuses() const {
    return std::make_pair(focus_1_, focus_2_);
  }

  double eccentricity() const {
    return c_ / a_;
  }

  Point center() const {
    Vec focuses = Vec(focus_1_, focus_2_);
    focuses /= 2;
    return focus_1_ + focuses;
  }

  std::pair<Line, Line> directrices() const {
    Point center = this->center();
    Vec vec = Vec(focus_1_, focus_2_);
    vec /= 2 * c_;
    vec *= (a_ * a_) / c_;
    Point d1 = center + vec;
    Point d2 = center - vec;
    Vec guide = vec.ort();
    Line l1 = Line(d1, guide);
    Line l2 = Line(d2, guide);
    return std::make_pair(l1, l2);
  }

  double perimeter() const override {
    double answer = PI * (3 * (a_ + b_) - sqrt((3 * a_ + b_) * (a_ + 3 * b_)));
    return answer;
  }

  double area() const override {
    return PI * a_ * b_;
  }

  bool operator==(const Shape& another) const override {
    const auto* el = dynamic_cast<const Ellipse*>(&another);
    if (el == nullptr) {
      return false;
    }
    if (focus_1_ == el->focus_1_ && focus_2_ == el->focus_2_&&
        equal_double(a_, el->a_) && equal_double(b_, el->b_)) {
      return true;
    }
    if (focus_2_ == el->focus_1_ && focus_1_ == el->focus_2_&&
        equal_double(a_, el->a_) && equal_double(b_, el->b_)) {
      return true;
    }
    return false;
  }

  bool operator!=(const Shape& another) const override {
    return !(*this == another);
  }

  bool isCongruentTo(const Shape& another) const override {
    const auto* el = dynamic_cast<const Ellipse*>(&another);
    if (el == nullptr) {
      return false;
    }
    if (equal_double(a_, el->a_) && equal_double(b_, el->b_)) {
      return true;
    }
    return false;
  }

  bool isSimilarTo(const Shape& another) const override {
    const auto* el = dynamic_cast<const Ellipse*>(&another);
    if (el == nullptr) {
      return false;
    }
    if (equal_double(a_ * el->b_, b_ * el->a_)) {
      return true;
    }
    return false;
  }

  bool containsPoint(const Point& p) const override {
    Point c = center();
    double x_offset = p.x - c.x;
    x_offset *= x_offset;
    double y_offset = p.y - c.y;
    y_offset *= y_offset;
    if (x_offset / (a_ * a_) + y_offset / (b_ * b_) <= 1 + precision) {
      return true;
    }
    return false;
  }

  void rotate(const Point& center, double angle) override {
    focus_1_.rotate(center, angle);
    focus_2_.rotate(center, angle);
  }

  void reflect(const Point& center) override {
    p_reflect(focus_1_, center);
    p_reflect(focus_2_, center);
  }

  void reflect(const Line& axis) override {
    p_reflect(focus_1_, axis);
    p_reflect(focus_2_, axis);
  }

  void scale(const Point& center, double coefficient) override {
    p_scale(focus_1_, center, coefficient);
    p_scale(focus_2_, center, coefficient);
    a_ *= coefficient;
    b_ *= coefficient;
    c_ *= coefficient;
  }

 protected:
  Point focus_1_;
  Point focus_2_;
  double a_;
  double b_;
  double c_;
};

class Circle : public Ellipse {
 public:
  Circle(const Point& p, double radius) : Ellipse(p, p, 2.0 * radius) {}

  double radius() const {
    return a_;
  }
};

class Rectangle : public Polygon {
 public:
  Rectangle(const Point& p1, const Point& p2, double k) : Polygon() {
    vertices_.resize(4);
    vertices_[0] = p1;
    vertices_[2] = p2;
    if (k < 1) {
      k = 1 / k;
    }
    Vec v1 =  Vec(p1, p2);
    double hypo = v1.Length();
    double small_leg = hypo / std::sqrt(1 + k * k);
    double big_leg = k * small_leg;
    double x = (small_leg * small_leg) / hypo;
    v1 /= hypo;
    v1 *= x;
    double h = (small_leg * big_leg) / hypo;
    Vec ort = v1.ort();
    ort /= ort.Length();
    ort *= h;
    Point p0 = p1 + v1;
    Point vert = p0 + ort;
    if (SignVecProduct(Vec(p1, vert), v1)) {
      vertices_[1] = vert;
      p0  = p2 - v1;
      vertices_[3]  = p0 - ort;
    } else {
      ort *= -1;
      vertices_[1] = p0 + ort;
      p0 = p2 - v1;
      vertices_[3] = p0 - ort;
    }
  }

  Point center() {
    Vec vec(vertices_[0], vertices_[2]);
    vec /= 2;
    Point center = vertices_[0] + vec;
    return center;
  }

  std::pair<Line, Line> diagonals() {
    return std::make_pair(
        Line(vertices_[0], vertices_[2]),
        Line(vertices_[1], vertices_[3]));
  }
};

class Square : public Rectangle {
 public:
  Square(const Point& p1, const Point& p2) : Rectangle(p1, p2, 1) {};

  Circle circumscribedCircle() {
    Point center = this->center();
    double radius = Vec(center, vertices_[0]).Length();
    Circle answer(center, radius);
    return answer;
  }

  Circle inscribedCircle() {
    Point center = this->center();
    double radius = Vec(center, vertices_[0]).Length() / std::sqrt(2);
    Circle answer(center, radius);
    return answer;
  }
};

class Triangle : public Polygon {
 public:
  Triangle(const Point& p1, const Point& p2, const Point& p3) : Polygon(p1, p2, p3) {};

  Circle circumscribedCircle() {
    double x12 = vertices_[0].x - vertices_[1].x;
    double x23 = vertices_[1].x - vertices_[2].x;
    double x31 = vertices_[2].x - vertices_[0].x;

    double y12 = vertices_[0].y - vertices_[1].y;
    double y23 = vertices_[1].y - vertices_[2].y;
    double y31 = vertices_[2].y - vertices_[0].y;

    double z1 = vertices_[0].x * vertices_[0].x + vertices_[0].y * vertices_[0].y;
    double z2 = vertices_[1].x * vertices_[1].x + vertices_[1].y * vertices_[1].y;
    double z3 = vertices_[2].x * vertices_[2].x + vertices_[2].y * vertices_[2].y;

    double zx = y12 * z3 + y23 * z1 + y31 * z2;
    double zy = x12 * z3 + x23 * z1 + x31 * z2;
    double z = x12 * y31 - y12 * x31;

    Point center(-(zx / (2 * z)), zy / (2 * z));
    Vec radius(center, vertices_[0]);
    return {center, radius.Length()};
  }

  Circle inscribedCircle() {
    Vec v1(vertices_[0], vertices_[1]);
    Vec v2(vertices_[0], vertices_[2]);

    Vec v3(vertices_[2], vertices_[0]);
    Vec v4(vertices_[2], vertices_[1]);

    v1.Normalize();
    v2.Normalize();
    v3.Normalize();
    v4.Normalize();

    Vec bis1vec = v1 + v2;
    Vec bis2vec = v3 + v4;

    Line bis1line(vertices_[0], bis1vec);
    Line bis2line(vertices_[2], bis2vec);

    Point center(Intersection(bis1line, bis2line));

    Line segment(vertices_[0], vertices_[2]);
    double radius = DistanceToLine(center, segment);
    return {center, radius};
  }

  Point centroid() {
    Vec v1(vertices_[1], vertices_[2]);
    v1 /= 2;
    Point v1center = vertices_[1] + v1;

    Line med1line(vertices_[0], v1center);

    Vec v2(vertices_[0], vertices_[1]);
    v2 /= 2;
    Point v2center = vertices_[0] + v2;

    Line med2line(vertices_[2], v2center);

    return Intersection(med1line, med2line);
  }

  Point orthocenter() {
    Vec v1(vertices_[0], vertices_[2]);
    Vec v1ort = v1.ort();

    Line h1(vertices_[1], v1ort);

    Vec v2(vertices_[0], vertices_[1]);
    Vec v2ort = v2.ort();

    Line h2(vertices_[2], v2ort);

    return Intersection(h1, h2);
  }

  Line EulerLine() {
    Point ort = orthocenter();
    Circle circle = circumscribedCircle();

    return {ort, circle.center()};
  }

  Circle ninePointsCircle() {
    Circle circumscribed_circle = circumscribedCircle();
    Point ort = orthocenter();

    Vec vec(ort, circumscribed_circle.center());
    vec /= 2;
    Point center = ort + vec;

    return {center, circumscribed_circle.radius() / 2};
  }
};

