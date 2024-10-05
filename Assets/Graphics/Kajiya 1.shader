Shader "Custom/KajiyaWithShadow"
{
    Properties
    {
        _DiffuseTerm("DiffuseTerm", Range(0.0, 1.0)) = 1
        _SpecularTerm("SpecularTerm", Range(0.0, 1.0)) = 1
        _SpecularPower("SpecularPower", Range(0.0, 1000.0)) = 10
        _Color ("Color", Color) = (1,1,1,1)
        _MainTex ("Texture", 2D) = "white" {}
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        //LOD 200
        
        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #include "UnityCG.cginc"
            #include "AutoLight.cginc"
            
            float _DiffuseTerm;
            float _SpecularTerm;
            float _SpecularPower;
            fixed4 _Color;
            sampler2D _MainTex;
            
            struct appdata
            {
                float4 vertex : POSITION;
                float3 tangent : TEXCOORD1;
            };
            
            struct v2f
            {
                //float2 uv : TEXCOORD0;
                float4 pos : SV_POSITION;
                float3 tangent : TANGENT;
                float4 lightDir : TEXCOORD1;
                float4 cameraDir : TEXCOORD2;
            };
            
            v2f vert (appdata v)
            {
                v2f o;
                o.pos = UnityObjectToClipPos(v.vertex);
                o.tangent = normalize(v.tangent);
                o.lightDir = mul(unity_WorldToObject, _WorldSpaceLightPos0);
                o.cameraDir = normalize(
                    mul(unity_WorldToObject, _WorldSpaceCameraPos)
                    - v.vertex); 
                return o;
            }
            
            float4 frag (v2f i) : SV_Target
            {
                float3 diffuse = _Color * _DiffuseTerm;

                float3 halfVec = normalize(i.lightDir + i.cameraDir);
                
                float3 specular = pow( 
                    sqrt(1 - dot(i.tangent, halfVec) * dot(i.tangent, halfVec)),
                    _SpecularPower) * _SpecularTerm;
                  
                float3 color = diffuse + specular;
                
                UNITY_LIGHT_ATTENUATION(atten, i, i.pos.xyz);

                //return float4(atten, 0, 0, 1);
                return float4(color, 1);
            }
            ENDCG
        }
        /*
        //Pass
        //{
            // Surface Shader 사용
            Tags { "LightMode" = "ForwardBase" }
            CGPROGRAM
            #pragma surface surf Standard fullforwardshadows

            sampler2D _MainTex;

            struct Input
            {
                float2 uv_MainTex;
            };

            // Surface 함수: albedo와 normal을 정의
            void surf(Input IN, inout SurfaceOutputStandard o)
            {
                half4 texColor = tex2D(_MainTex, IN.uv_MainTex);
                o.Albedo = texColor.rgb;  // 기본 색상
                o.Alpha = texColor.a;     // 투명도 설정
            }
            ENDCG
        //}
        */
    }
    FallBack "Diffuse"
}
