
# Copyright (c) 2021-2024, University of Luxembourg
# Copyright (c) 2021-2024, Heinrich-Heine University Duesseldorf
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

@testset "Automated QUality Assurance" begin
    # We can't do Aqua.test_all here (yet) because the ambiguity tests fail in
    # deps. Instead let's pick the other sensible tests.
    Aqua.test_deps_compat(COBREXA)
    Aqua.test_project_extras(COBREXA)
    Aqua.test_project_toml_formatting(COBREXA)
    #Aqua.test_stale_deps(COBREXA) # currently seems broken
    Aqua.test_unbound_args(COBREXA)
    Aqua.test_undefined_exports(COBREXA)
end
